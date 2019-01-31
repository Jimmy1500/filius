#include "Swaption.h"

Swaption::Swaption()
    : RateInstrument(RateInstrumentType::RIT_SWAPTION, "Swaption with "),
    Model(nullptr),
    Params(new size_t[SWPT::NUM_PARAMS]()),
    Flags (new size_t[SWPT::NUM_PARAMS]),
    Dirty(0),
    Prices(nullptr),
    Payoff(0.)
{
    Flags[SWPT::NTL] = SWPT::GET_VALUE;
    Flags[SWPT::STK] = SWPT::GET_PAYOFF;
    Flags[SWPT::STL] = SWPT::GET_ZCBP;
}

Swaption::Swaption(RateModel* model)
    : RateInstrument(RateInstrumentType::RIT_SWAPTION, "Swaption with "),
    Model(model),
    Params(new size_t[SWPT::NUM_PARAMS]()),
    Flags (new size_t[SWPT::NUM_PARAMS]),
    Dirty(0),
    Prices(nullptr),
    Payoff(0.)
{
    if (Model){ Description.append( Model->getModelDescription()); }
    Flags[SWPT::NTL] = SWPT::GET_VALUE;
    Flags[SWPT::STK] = SWPT::GET_PAYOFF;
    Flags[SWPT::STL] = SWPT::GET_ZCBP;
}

Swaption::Swaption(Swaption & swaption)
    : RateInstrument(RateInstrumentType::RIT_SWAPTION, "Swaption with "),
    Model(swaption.getModel()),
    Params(new size_t[SWPT::NUM_PARAMS]()),
    Flags (new size_t[SWPT::NUM_PARAMS]),
    Dirty(0),
    Prices(nullptr),
    Payoff(swaption.Payoff)
{
    if (Model){ Description.append( Model->getModelDescription()); }
    size_t i;
    for ( i = 0; i < SWPT::NUM_PARAMS; ++i ){
        Params[i] = swaption.Params[i];
        Flags[i] = swaption.Flags[i];
    }
}

Swaption::~Swaption(){
    delete [] Flags;
    delete [] Params;
    delete Prices;
}

double Swaption::getModelValue(){
    if (!Model){
#ifdef __DEBUG__
        DEBUG("Swaption: no interest model detected, defaulting to presumed market value")
#endif
        return MarketValue;
    }
    switch (Model->getModelType()){
        case RMT_G2PP:
            {
                if ( G2PP * g2pp = dynamic_cast<G2PP *>(Model) ){
                    size_t  nthreads = g2pp->getPeriphery(G2::NTHREADS);
                    size_t  npaths   = g2pp->getPeriphery(G2::NPATHS);
                    double  strike   = Params[SWPT::STK];
                    double  settle   = Params[SWPT::STL];
                    double *terms    = Terms;
                    size_t  nterms   = NumTerms;
#ifdef __REGEN__
                    markDirtyAll();
#else
                    if ( g2pp->isDirty(G2::GENERATION) || g2pp->isDirty(G2::EVOLUTION) ) {
                        markDirtyFrom(SWPT::GET_ZCBP);
                    }
#endif

                    if (Prices){
                        if (Prices->resize(nterms, npaths)){ markDirtyAll(); }
                    }else{
                        Prices = new mat2d(nterms, npaths); markDirtyAll();
                    }

                    if (isDirty(SWPT::GET_ZCBP)){
                        g2pp->getZCBP(Prices->value, settle, terms, nterms, npaths);
                        clearDirty(SWPT::GET_ZCBP);
                    }

                    if (isDirty(SWPT::GET_PAYOFF)){
                        Payoff = 0.0;
                        switch(nthreads){
                            case 1:
                                {
                                    size_t path, term;
                                    double float_leg, fix_leg, fix_factor;
                                    for (path = 0 ; path < Prices->npaths; ++path){
                                        fix_factor = ( terms[0] - settle ) * Prices->value[path][0];
                                        for (term = 1; term < Prices->nterms; ++term){
                                            fix_factor += ( terms[term] - terms[term-1] ) * Prices->value[path][term];
                                        }
                                        fix_leg = strike * fix_factor;
                                        float_leg = Prices->value[path][0] - Prices->value[path][nterms-1];
                                        Payoff += max(float_leg - fix_leg);
                                    }
                                    break;
                                }// case 1
                            default:
                                {
                                    size_t  i;
                                    size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads;
                                    thread * threads = new thread[nthreads];
                                    for (i = 0; i < last_thread; ++i){
                                        threads[i] = thread([this,i,paths_per_thread,strike,settle,terms,nterms]()->void{
                                            size_t begin = i*paths_per_thread;
                                            size_t end   = (i+1)*paths_per_thread;
                                            size_t path, term;
                                            double float_leg, fix_leg, fix_factor, pay_off = 0.;
                                            for (path = begin; path < end; ++path){
                                                fix_factor = ( terms[0] - settle ) * Prices->value[path][0];
                                                for (term = 1; term < Prices->nterms; ++term){
                                                    fix_factor += ( terms[term] - terms[term-1] ) * Prices->value[path][term];
                                                }
                                                fix_leg = strike * fix_factor;
                                                float_leg = Prices->value[path][0] - Prices->value[path][nterms-1];
                                                pay_off += max(float_leg - fix_leg);
                                            }
                                            mtx.lock(); Payoff += pay_off; mtx.unlock();
                                        });
                                    }

                                    threads[last_thread] = thread([this,npaths,last_thread,paths_per_thread,strike,settle,terms,nterms]()->void{
                                        size_t begin = last_thread*paths_per_thread;
                                        size_t path, term;
                                        double float_leg, fix_leg, fix_factor, pay_off = 0.;
                                        for (path = begin; path < npaths; ++path){
                                            fix_factor = ( terms[0] - settle ) * Prices->value[path][0];
                                            for (term = 1; term < Prices->nterms; ++term){
                                                fix_factor += ( terms[term] - terms[term-1] ) * Prices->value[path][term];
                                            }
                                            fix_leg = strike * fix_factor;
                                            float_leg = Prices->value[path][0] - Prices->value[path][nterms-1];
                                            pay_off += max(float_leg - fix_leg);
                                        }
                                        mtx.lock(); Payoff += pay_off; mtx.unlock();
                                    });

                                    for (i = 0; i < nthreads; ++i){ threads[i].join(); }

                                    thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                                    garbage_collection.detach();
                                    break;
                                }// default
                        }// switch (nthreads)
                        clearDirty(SWPT::GET_PAYOFF);
                    }// if (isDirty(SWPT::GET_PAYOFF))

                    if (isDirty(SWPT::GET_VALUE)){
                        Payoff *= Params[SWPT::NTL] / (double)npaths;
                        clearDirty(SWPT::GET_VALUE);
                    }

                    return Payoff;

                } else {// if ( G2PP * g2pp = dynamic_cast<G2PP *>(Model) )
                    REPORT_ERROR(Model->ModelError, "Unable to resolve rate model type (G2PP)", 0)
                }

                break;
            }
        case RMT_BLACK:
            {
                return MarketValue;
                break;
            }
        default:
            {
                return MarketValue;
                break;
            }
    }
}

double Swaption::getModelValue(double notional, double strike, double settlement, double * terms, size_t nterms){
    size_t keys[SWPT::NUM_PARAMS] = {SWPT::NTL, SWPT::STK, SWPT::STL};
    double values[SWPT::NUM_PARAMS] = {notional, strike, settlement};
    setParameters(keys, values, SWPT::NUM_PARAMS);
    setTerms(terms, nterms);
    return getModelValue();
}
