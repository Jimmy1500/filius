#include "G2PP.h"

G2PP::G2PP()
    : RateModel(RateModelType::RMT_G2PP, "Gaussian 2-factor short rate model(G2++)"),
    Sim(new Simulation()),
    YieldCurve(nullptr),
    Samples(nullptr),
    Factors(nullptr),
    Results(nullptr),
    Settlement(0.0),
    Maturity(0.0),
    ZCB(0.0),
    Dirty(0),
    Flags(new size_t[G2::NUM_PARAMS]),
    Coefs(new double[G2::NUM_PARAMS]()),
    Peris(new size_t[G2::NUM_PERIS]())
{
    fill(Flags, Flags+G2::NUM_PARAMS, G2::EVOLUTION);
    this->markDirtyAll();
}

G2PP::G2PP(Curve * curve)
    : RateModel(RateModelType::RMT_G2PP, "G2++: Gaussian 2-factor short rate model"),
    Sim(new Simulation()),
    YieldCurve(curve),
    Samples(nullptr),
    Factors(nullptr),
    Results(nullptr),
    Settlement(0.0),
    Maturity(0.0),
    ZCB(0.0),
    Dirty(0),
    Flags(new size_t[G2::NUM_PARAMS]),
    Coefs(new double[G2::NUM_PARAMS]()),
    Peris(new size_t[G2::NUM_PERIS]())
{
    fill(Flags, Flags+G2::NUM_PARAMS, G2::EVOLUTION);
    this->markDirtyAll();
}

G2PP::~G2PP(){
    delete [] Flags;
    delete [] Coefs;
    delete [] Peris;
    delete Sim;
    delete Samples;
    delete Factors;
    delete Results;
    YieldCurve = nullptr;//YieldCurve won't be deleted to allow external persistence pattern, and centralized instantiation and destruction pattern
}

//----Calculation Services-----
double G2PP::getZCBP(double t, double T){
    if (!YieldCurve){ REPORT_ERROR(ModelError, "Current market condition not set", 0) }
    if ( T < t ){ REPORT_ERROR(ModelError, "Maturity is prior to settlement", 0) }
    else if ( t == T ){ return 1.0; }
    else if ( t == 0.0 ){ return YieldCurve->P(T); }
    else{
        if ( Settlement != t ){ Settlement = t; markDirtyFrom(G2::EVOLUTION); }
        if ( Maturity != T ){ Maturity = T; markDirtyFrom(G2::SIMULATION); }
    }
    
    size_t npaths = Peris[G2::NPATHS];
    size_t ndims = Peris[G2::NDIMS];
    size_t nthreads = Peris[G2::NTHREADS];
    if (!npaths || !ndims){ REPORT_ERROR(ModelError, "Simulation periphery not set", 0) }
    
    if (nthreads < 1){
        nthreads = thread::hardware_concurrency();
#ifdef __DEBUG__
        WARN("Simulation threads number illogical, default to cpu core population")
#endif
    }
    if (nthreads > npaths){
        nthreads = npaths;
#ifdef __DEBUG__
        WARN("Simulation threads number unreasonable, default to simulation periphery(npaths)")
#endif
    }
    
    if (Samples){
        if (Samples->resize(1,npaths,ndims)){ markDirtyAll(); }
    }else{
        Samples = new mat3d(1,npaths,ndims); markDirtyAll();
    }
    if (Factors){
        if (Factors->resize(npaths,ndims)){ markDirtyFrom(G2::EVOLUTION); }
    }else{
        Factors = new mat2d(npaths,ndims); markDirtyFrom(G2::EVOLUTION);
    }

    if (isDirty(G2::GENERATION)){ Sim->ngen(*Samples,nthreads); clearDirty(G2::GENERATION); }
    
    double a = Coefs[G2::A];
    double b = Coefs[G2::B];
    double vol1 = Coefs[G2::SIGMA_1];
    double vol2 = Coefs[G2::SIGMA_2];
    double rho = Coefs[G2::RHO];
    double T_t = T-t;

    double * X = Factors->value[G2::X];
    double * Y = Factors->value[G2::Y];

    size_t i;
    if (isDirty(G2::EVOLUTION)){
        if (Results){ Results->resize(nthreads); markDirtyFrom(G2::SIMULATION); }
        else{ Results = new mat1d(nthreads); markDirtyFrom(G2::SIMULATION); }

        double ** Z_X = Samples->value[G2::X];
        double ** Z_Y = Samples->value[G2::Y];

        double chol_y = sqrt(1.-rho*rho);

        switch (nthreads){
            case 1:
                {
                    size_t path;
                    for (path = 0; path < npaths; ++path){
                        if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                        else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                        if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                        else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                        Results->value[0] += .5*( exp(-M_XY( X[path], Y[path], T_t)) + exp(-M_XY(-X[path],-Y[path], T_t)) );
                    }
                    break;
                }
            default:
                {
                    size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads;

                    thread * threads = new thread[nthreads];
                    for (i = 0; i < last_thread; ++i){
                        threads[i] = thread([this,i,paths_per_thread,a,b,vol1,vol2,rho,T_t,chol_y,t,Z_X,Z_Y,X,Y]()->void{
                            size_t begin = i*paths_per_thread;
                            size_t end   = (i+1)*paths_per_thread;
                            size_t path;
                            for (path = begin; path < end; ++path){
                                if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                                else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                                if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                                else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                                Results->value[i] += .5*( exp(-M_XY( X[path], Y[path], T_t)) + exp(-M_XY(-X[path],-Y[path], T_t)) );
                            }
                        });
                    }

                    threads[last_thread] = thread([this,npaths,last_thread,paths_per_thread,a,b,vol1,vol2,rho,T_t,chol_y,t,Z_X,Z_Y,X,Y]()->void{
                        size_t begin = last_thread*paths_per_thread;
                        size_t path;
                        for (path = begin; path < npaths; ++path){
                            if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                            else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                            if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                            else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                            Results->value[last_thread] += .5*( exp(-M_XY( X[path], Y[path], T_t)) + exp(-M_XY(-X[path],-Y[path], T_t)) );
                        }
                    });

                    for (i = 0; i < nthreads; ++i){
                        threads[i].join();
                    }

                    thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                    garbage_collection.detach();
                    break;
                }//default
        }//switch(threads)
        clearDirty(G2::EVOLUTION);
        clearDirty(G2::SIMULATION); //isDirty(G2::EVOLUTION) => isDirty(G2::SIMULATION), but NOT vise versa
    }

    if (isDirty(G2::SIMULATION)){
        if (Results){ Results->resize(nthreads); markDirtyFrom(G2::SIMULATION);
        }else{ Results = new mat1d(nthreads); markDirtyFrom(G2::SIMULATION); }

        size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads;
        thread * threads = new thread[nthreads];

        for (i = 0; i < last_thread; ++i){
            threads[i] = thread([this,i,paths_per_thread,a,b,T_t,X,Y]()->void{
                size_t begin = i*paths_per_thread;
                size_t end   = (i+1)*paths_per_thread;
                size_t path;
                for (path = begin; path < end; ++path){
                    Results->value[i] += .5*( exp(-M_XY( X[path], Y[path], T_t)) + exp(-M_XY(-X[path],-Y[path], T_t)) );
                }
            });
        }

        threads[last_thread] = thread([this,npaths,last_thread,paths_per_thread,a,b,T_t,X,Y]()->void{
            size_t begin = last_thread*paths_per_thread;
            size_t path;
            for (path = begin; path < npaths; ++path){
                Results->value[last_thread] += .5*( exp(-M_XY( X[path], Y[path], T_t)) + exp(-M_XY(-X[path],-Y[path], T_t)) );
            }
        });

        for (i = 0; i < nthreads; ++i){
            threads[i].join();
        }

        thread garbage_collection = thread([threads]()->void{ delete [] threads; });
        garbage_collection.detach();

        clearDirty(G2::SIMULATION);
    }

    if (isDirty(G2::CALCULATION)){
        ZCB = 0.0;
        for (i = 0; i < nthreads; ++i){
            ZCB += Results->value[i];
        }

        ZCB *= YieldCurve->P(t,T) * exp( .5*(V_XY(T_t)-V_XY(T)+V_XY(t)) ) / (double)npaths;
        clearDirty(G2::CALCULATION);
    }

    return ZCB;
}

void G2PP::getZCBP(double * prices, double * t, double T, size_t nterms){
    if ( !prices ){ REPORT_ERROR(ModelError, "Zero coupon bond output container corrupted", 0) }
    else if ( !nterms ){ REPORT_ERROR(ModelError, "User indicates output container size is zero, calculation will not commence", 0) }
    else if ( prices[0] != 0.0 || prices[nterms-1] != 0.0 ){ fill(prices, prices+nterms, 0.0); }
    else{
        size_t npaths = Peris[G2::NPATHS];
        size_t ndims = Peris[G2::NDIMS];
        size_t nthreads = Peris[G2::NTHREADS];

        if (nthreads < 1){
            nthreads = thread::hardware_concurrency();
#ifdef __DEBUG__
            WARN("Simulation threads number illogical, default to cpu core population")
#endif
        }
        if (nthreads > npaths){
            nthreads = npaths;
#ifdef __DEBUG__
            WARN("Simulation threads number unreasonable, default to simulation periphery(npaths)")
#endif
        }

#ifdef __REGEN__
        if (Samples){
            if (Samples->resize(nterms,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(nterms,npaths,ndims); markDirtyAll();
        }
#else
        if (Samples){
            if (Samples->resize(1,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(1,npaths,ndims); markDirtyAll();
        }
#endif
        if (isDirty(G2::GENERATION)){ Sim->ngen(*Samples,nthreads); clearDirty(G2::GENERATION); }
        double a = Coefs[G2::A];
        double b = Coefs[G2::B];
        double vol1 = Coefs[G2::SIGMA_1];
        double vol2 = Coefs[G2::SIGMA_2];
        double rho = Coefs[G2::RHO];
        double chol_y = sqrt(1.-rho*rho);

        switch (nthreads){
            case 1:
                {
                    double t_,T_t,Z_X,Z_Y,X,Y;
                    size_t term, path;
                    for (term = 0; term < nterms; ++term){
                        t_ = t[term]; T_t = T-t_;
                        for (path = 0; path < npaths; ++path){
#ifdef __REGEN__
                            Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                            Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                            if (a==0.0){ X = vol1*sqrt(t_)*Z_X; }
                            else{        X = vol1*sqrt(.5*(1.-exp(-2.*a*t_))/a)*Z_X; }
                            if (b==0.0){ Y = vol2*sqrt(t_)*(rho*Z_X-chol_y*Z_Y); }
                            else{        Y = vol2*sqrt(.5*(1.-exp(-2.*b*t_))/b)*(rho*Z_X-chol_y*Z_Y); }
                            prices[term] += .5*( exp(-M_XY( X, Y, T_t)) + exp(-M_XY(-X,-Y, T_t)) );
                        }
                        prices[term] *= YieldCurve->P(t_,T) * exp( .5*(V_XY(T_t)-V_XY(T)+V_XY(t_)) ) / (double)npaths;
                    }
                    break;
                }
            default:
                {
                    size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads, i;
                    thread * threads = new thread[nthreads];
                    for (i = 0; i < last_thread; ++i){
                        threads[i] = thread([this,nterms,i,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,prices]()->void{
                            double t_,T_t,Z_X,Z_Y,X,Y,coef,mean;
                            size_t term, path, begin, end;
                            for (term = 0; term < nterms; ++term){
                                t_ = t[term]; T_t = T-t_; coef = YieldCurve->P(t_,T) * exp( .5*(V_XY(T_t)-V_XY(T)+V_XY(t_)) ) / (double)npaths; mean = 0.0;
                                begin = i*paths_per_thread, end = (i+1)*paths_per_thread;
                                for (path = begin; path < end; ++path){
#ifdef __REGEN__
                                    Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                    Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                    if (a==0){ X = vol1*sqrt(t_)*Z_X; }
                                    else{      X = vol1*sqrt(.5*(1.-exp(-2.*a*t_))/a)*Z_X; }
                                    if (b==0){ Y = vol2*sqrt(t_)*(rho*Z_X-chol_y*Z_Y); }
                                    else{      Y = vol2*sqrt(.5*(1.-exp(-2.*b*t_))/b)*(rho*Z_X-chol_y*Z_Y); }
                                    mean += coef*.5*( exp(-M_XY( X, Y, T_t)) + exp(-M_XY(-X,-Y, T_t)) );
                                }
                                mtx.lock(); prices[term] += mean; mtx.unlock();
                            }
                        });
                    }

                    threads[last_thread] = thread([this,nterms,last_thread,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,prices]()->void{
                        double t_,T_t,Z_X,Z_Y,X,Y,coef,mean;
                        size_t term, path, begin;
                        for (term = 0; term < nterms; ++term){
                            t_ = t[term]; T_t = T-t_; coef = YieldCurve->P(t_,T) * exp( .5*(V_XY(T_t)-V_XY(T)+V_XY(t_)) ) / (double)npaths; mean = 0.0;
                            begin = last_thread*paths_per_thread;
                            for (path = begin; path < npaths; ++path){
#ifdef __REGEN__
                                Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                if (a==0){ X = vol1*sqrt(t_)*Z_X; }
                                else{      X = vol1*sqrt(.5*(1.-exp(-2.*a*t_))/a)*Z_X; }
                                if (b==0){ Y = vol2*sqrt(t_)*(rho*Z_X-chol_y*Z_Y); }
                                else{      Y = vol2*sqrt(.5*(1.-exp(-2.*b*t_))/b)*(rho*Z_X-chol_y*Z_Y); }
                                mean += coef*.5*( exp(-M_XY( X, Y, T_t)) + exp(-M_XY(-X,-Y, T_t)) );
                            }
                            mtx.lock(); prices[term] += mean; mtx.unlock();
                        }
                    });

                    for (i = 0; i < nthreads; ++i){
                        threads[i].join();
                    }

                    thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                    garbage_collection.detach();
                    break;
                }//default (multithreaded)
        }//switch(nthreads)
    }//if(nterms>0)
}

void G2PP::getZCBP(double * prices, double t, double * T, size_t nterms){
    if ( !prices ){ REPORT_ERROR(ModelError, "Zero coupon bond output container corrupted", 0) }
    else if ( !nterms ){ REPORT_ERROR(ModelError, "User indicates output container size is zero, calculation will not commence", 0) }
    else if ( prices[0] != 0.0 || prices[nterms-1] != 0.0 ){ fill(prices, prices+nterms, 0.0); }
    if ( Settlement != t ){ Settlement = t; markDirtyFrom(G2::EVOLUTION); }
    else{
        size_t ndims    = Peris[G2::NDIMS];
        size_t npaths   = Peris[G2::NPATHS];
        size_t nthreads = Peris[G2::NTHREADS];

        if (nthreads < 1){
            nthreads = thread::hardware_concurrency();
#ifdef __DEBUG__
            WARN("Simulation threads number illogical, default to cpu core population")
#endif
        }
        if (nthreads > npaths){
            nthreads = npaths;
#ifdef __DEBUG__
            WARN("Simulation threads number unreasonable, default to simulation periphery(npaths)")
#endif
        }
#ifdef __REGEN__
        if (Samples){
            if (Samples->resize(nterms,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(nterms,npaths,ndims); markDirtyAll();
        }
#else
        if (Samples){
            if (Samples->resize(1,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(1,npaths,ndims); markDirtyAll();
        }
#endif
        if (Factors){
            if (Factors->resize(npaths,ndims)){ markDirtyFrom(G2::EVOLUTION); }
        }else{
            Factors = new mat2d(npaths,ndims); markDirtyFrom(G2::EVOLUTION);
        }

        if (isDirty(G2::GENERATION)){ Sim->ngen(*Samples,nthreads); clearDirty(G2::GENERATION); }
        
        if (isDirty(G2::EVOLUTION)){
            double a = Coefs[G2::A];
            double b = Coefs[G2::B];
            double vol1 = Coefs[G2::SIGMA_1];
            double vol2 = Coefs[G2::SIGMA_2];
            double rho = Coefs[G2::RHO];
            double chol_y = sqrt(1.-rho*rho);

            double * X = Factors->value[G2::X];
            double * Y = Factors->value[G2::Y];
            switch (nthreads){
                case 1:
                    {
                        double T_,T_t,Z_X,Z_Y,coef;
                        size_t term, path;
                        for (term = 0; term < nterms; ++term){
                            T_= T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths;
                            for (path = 0; path < npaths; ++path){
#ifdef __REGEN__
                                Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                if (a==0.0){ X[path] = vol1*sqrt(t)*Z_X; }
                                else{        X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                if (b==0.0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                else{        Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                prices[term] += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                            }
                        }
                        break;
                    }
                default:
                    {
                        size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads, i;
                        thread * threads = new thread[nthreads];
                        for (i = 0; i < last_thread; ++i){
                            threads[i] = thread([this,nterms,i,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                                double T_,T_t,Z_X,Z_Y,coef,mean;
                                size_t term, path, begin = i*paths_per_thread, end = (i+1)*paths_per_thread;
                                for (term = 0; term < nterms; ++term){
                                    T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths; mean = 0;
                                    for (path = begin; path < end; ++path){
#ifdef __REGEN__
                                        Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                        Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                        if (a==0){ X[path] = vol1*sqrt(t)*Z_X; }
                                        else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                        if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                        else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                        mean += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                    }
                                    mtx.lock(); prices[term] += mean; mtx.unlock();
                                }
                            });
                        }

                        threads[last_thread] = thread([this,nterms,last_thread,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                            double T_,T_t,Z_X,Z_Y,coef,mean;
                            size_t term, path, begin = last_thread*paths_per_thread;
                            for (term = 0; term < nterms; ++term){
                                T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths; mean = 0;
                                for (path = begin; path < npaths; ++path){
#ifdef __REGEN__
                                    Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                    Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                    if (a==0){ X[path] = vol1*sqrt(t)*Z_X; }
                                    else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                    if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                    else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                    mean += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                }
                                mtx.lock(); prices[term] += mean; mtx.unlock();
                            }
                        });

                        for (i = 0; i < nthreads; ++i){
                            threads[i].join();
                        }

                        thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                        garbage_collection.detach();
                        break;
                    }// default
            }// switch (nthreads)
            clearDirty(G2::EVOLUTION);
        } else {// if (isDirty(G2::EVOLUTION))
            double a = Coefs[G2::A];
            double b = Coefs[G2::B];
            double vol1 = Coefs[G2::SIGMA_1];
            double vol2 = Coefs[G2::SIGMA_2];
            double rho = Coefs[G2::RHO];
            double chol_y = sqrt(1.-rho*rho);

            double * X = Factors->value[G2::X];
            double * Y = Factors->value[G2::Y];
            switch (nthreads){
                case 1:
                    {
                        double T_,T_t,coef;
                        size_t term, path;
                        for (term = 0; term < nterms; ++term){
                            T_= T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths;
                            for (path = 0; path < npaths; ++path){
                                prices[term] += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                            }
                        }
                        break;
                    }
                default:
                    {
                        size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads, i;
                        thread * threads = new thread[nthreads];
                        for (i = 0; i < last_thread; ++i){
                            threads[i] = thread([this,nterms,i,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                                double T_,T_t,coef,mean;
                                size_t term, path, begin = i*paths_per_thread, end = (i+1)*paths_per_thread;
                                for (term = 0; term < nterms; ++term){
                                    T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths; mean = 0;
                                    for (path = begin; path < end; ++path){
                                        mean += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                    }
                                    mtx.lock(); prices[term] += mean; mtx.unlock();
                                }
                            });
                        }

                        threads[last_thread] = thread([this,nterms,last_thread,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                            double T_,T_t,coef,mean;
                            size_t term, path, begin = last_thread*paths_per_thread;
                            for (term = 0; term < nterms; ++term){
                                T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) ) / (double)npaths; mean = 0;
                                for (path = begin; path < npaths; ++path){
                                    mean += coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                }
                                mtx.lock(); prices[term] += mean; mtx.unlock();
                            }
                        });

                        for (i = 0; i < nthreads; ++i){
                            threads[i].join();
                        }

                        thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                        garbage_collection.detach();
                        break;
                    }// default
            }// switch (nthreads)
        }// if (!isDirty(G2::EVOLUTION))
    }// if (nterms)
}

void G2PP::getZCBP(double ** prices, double t, double * T, size_t nterms, size_t npaths){
    if ( !prices ){ REPORT_ERROR(ModelError, "Zero coupon bond output container corrupted", 0) }
    else if ( !nterms || !npaths ){ REPORT_ERROR(ModelError, "User indicates output container size is zero, calculation will not commence", 0) }
    else if ( npaths != Peris[G2::NPATHS] ){ REPORT_ERROR(ModelError, "Simulation output periphery doesn't match user specification: npaths", 0) }
    else {
        if ( Settlement != t ){ Settlement = t; markDirtyFrom(G2::EVOLUTION); }

        size_t i;
        for ( i = 0; i < npaths; ++i ){ 
            if ( prices[i][0] != 0.0 || prices[i][nterms-1] != 0.0 ){
                fill(prices[i], prices[i]+nterms, 0.0); 
            }
        }

        size_t ndims    = Peris[G2::NDIMS];
        size_t nthreads = Peris[G2::NTHREADS];

        if (nthreads < 1){
            nthreads = thread::hardware_concurrency();
#ifdef __DEBUG__
            WARN("Simulation threads number illogical, default to cpu core population")
#endif
        }
        if (nthreads > npaths){
            nthreads = npaths;
#ifdef __DEBUG__
            WARN("Simulation threads number unreasonable, default to simulation periphery(npaths)")
#endif
        }
#ifdef __REGEN__
        if (Samples){
            if (Samples->resize(nterms,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(nterms,npaths,ndims); markDirtyAll();
        }
#else
        if (Samples){
            if (Samples->resize(1,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(1,npaths,ndims); markDirtyAll();
        }
#endif
        if (Factors){
            if (Factors->resize(npaths,ndims)){ markDirtyFrom(G2::EVOLUTION); }
        }else{
            Factors = new mat2d(npaths,ndims); markDirtyFrom(G2::EVOLUTION);
        }

        if (isDirty(G2::GENERATION)){ Sim->ngen(*Samples,nthreads); clearDirty(G2::GENERATION); }
        
        if (isDirty(G2::EVOLUTION)){
            double a = Coefs[G2::A];
            double b = Coefs[G2::B];
            double vol1 = Coefs[G2::SIGMA_1];
            double vol2 = Coefs[G2::SIGMA_2];
            double rho = Coefs[G2::RHO];
            double chol_y = sqrt(1.-rho*rho);

            double * X = Factors->value[G2::X];
            double * Y = Factors->value[G2::Y];
            switch (nthreads){
                case 1:
                    {
                        double T_,T_t,Z_X,Z_Y,coef;
                        size_t term, path;
                        for (term = 0; term < nterms; ++term){
                            T_= T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                            for (path = 0; path < npaths; ++path){
#ifdef __REGEN__
                                Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                if (a==0.0){ X[path] = vol1*sqrt(t)*Z_X; }
                                else{        X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                if (b==0.0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                else{        Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                            }
                        }
                        break;
                    }
                default:
                    {
                        size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads;
                        thread * threads = new thread[nthreads];
                        for (i = 0; i < last_thread; ++i){
                            threads[i] = thread([this,nterms,i,paths_per_thread,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                                double T_,T_t,Z_X,Z_Y,coef;
                                size_t term, path, begin = i*paths_per_thread, end = (i+1)*paths_per_thread;
                                for (term = 0; term < nterms; ++term){
                                    T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                                    for (path = begin; path < end; ++path){
#ifdef __REGEN__
                                        Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                        Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                        if (a==0){ X[path] = vol1*sqrt(t)*Z_X; }
                                        else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                        if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                        else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                        prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                    }
                                }
                            });
                        }

                        threads[last_thread] = thread([this,nterms,last_thread,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                            double T_,T_t,Z_X,Z_Y,coef;
                            size_t term, path, begin = last_thread*paths_per_thread;
                            for (term = 0; term < nterms; ++term){
                                T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                                for (path = begin; path < npaths; ++path){
#ifdef __REGEN__
                                    Z_X = Samples->value[G2::X][path][term]; Z_Y = Samples->value[G2::Y][path][term];
#else
                                    Z_X = Samples->value[G2::X][path][0]; Z_Y = Samples->value[G2::Y][path][0];
#endif
                                    if (a==0){ X[path] = vol1*sqrt(t)*Z_X; }
                                    else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X; }
                                    if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X-chol_y*Z_Y); }
                                    else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X-chol_y*Z_Y); }
                                    prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                }
                            }
                        });

                        for (i = 0; i < nthreads; ++i){
                            threads[i].join();
                        }

                        thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                        garbage_collection.detach();
                        break;
                    }// default
            }// switch (nthreads)
            clearDirty(G2::EVOLUTION);
        } else {// if (isDirty(G2::EVOLUTION))
            double a = Coefs[G2::A];
            double b = Coefs[G2::B];
            double vol1 = Coefs[G2::SIGMA_1];
            double vol2 = Coefs[G2::SIGMA_2];
            double rho = Coefs[G2::RHO];
            double chol_y = sqrt(1.-rho*rho);

            double * X = Factors->value[G2::X];
            double * Y = Factors->value[G2::Y];
            switch (nthreads){
                case 1:
                    {
                        double T_,T_t,coef;
                        size_t term, path;
                        for (term = 0; term < nterms; ++term){
                            T_= T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                            for (path = 0; path < npaths; ++path){
                                prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                            }
                        }
                        break;
                    }
                default:
                    {
                        size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads, i;
                        thread * threads = new thread[nthreads];
                        for (i = 0; i < last_thread; ++i){
                            threads[i] = thread([this,nterms,i,paths_per_thread,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                                double T_,T_t,coef;
                                size_t term, path, begin = i*paths_per_thread, end = (i+1)*paths_per_thread;
                                for (term = 0; term < nterms; ++term){
                                    T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                                    for (path = begin; path < end; ++path){
                                        prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                    }
                                }
                            });
                        }

                        threads[last_thread] = thread([this,nterms,last_thread,paths_per_thread,npaths,a,b,vol1,vol2,rho,chol_y,t,T,X,Y,prices]()->void{
                            double T_,T_t,coef;
                            size_t term, path, begin = last_thread*paths_per_thread;
                            for (term = 0; term < nterms; ++term){
                                T_ = T[term]; T_t = T_-t; coef = YieldCurve->P(t,T_) * exp( .5*(V_XY(T_t)-V_XY(T_)+V_XY(t)) );
                                for (path = begin; path < npaths; ++path){
                                    prices[path][term] = coef*.5*( exp(-M_XY(X[path], Y[path],T_t)) + exp(-M_XY(-X[path],-Y[path],T_t)) );
                                }
                            }
                        });

                        for (i = 0; i < nthreads; ++i){
                            threads[i].join();
                        }

                        thread garbage_collection = thread([threads]()->void{ delete [] threads; });
                        garbage_collection.detach();
                        break;
                    }// default
            }// switch (nthreads)
        }// if (!isDirty(G2::EVOLUTION))
    }// if (nterms)
}

mat2d * G2PP::getFactors(double t){
    if ( Settlement != t ){ Settlement = t; markDirtyFrom(G2::EVOLUTION); }
    if (isDirty(G2::EVOLUTION)){
        size_t npaths = Peris[G2::NPATHS];
        size_t ndims = Peris[G2::NDIMS];
        size_t nthreads = Peris[G2::NTHREADS];
        if (!npaths || !ndims){ REPORT_ERROR(ModelError, "Simulation periphery not set", 0) }
        if (nthreads < 1){
            nthreads = thread::hardware_concurrency();
#ifdef __DEBUG__
            WARN("Simulation threads number illogical, default to cpu core population")
#endif
        }
        if (nthreads > npaths){
            nthreads = npaths;
#ifdef __DEBUG__
            WARN("Simulation threads number unreasonable, default to simulation periphery(npaths)")
#endif
        }

        if (Samples){
            if (Samples->resize(1,npaths,ndims)){ markDirtyAll(); }
        }else{
            Samples = new mat3d(1,npaths,ndims); markDirtyAll();
        }
        if (Factors){
            if (Factors->resize(npaths,ndims)){ markDirtyFrom(G2::EVOLUTION); }
        }else{
            Factors = new mat2d(npaths,ndims); markDirtyFrom(G2::EVOLUTION);
        }

        if (isDirty(G2::GENERATION)){ Sim->ngen(*Samples,nthreads); clearDirty(G2::GENERATION); }

        double a = Coefs[G2::A];
        double b = Coefs[G2::B];
        double vol1 = Coefs[G2::SIGMA_1];
        double vol2 = Coefs[G2::SIGMA_2];
        double rho = Coefs[G2::RHO];

        double * X = Factors->value[G2::X];
        double * Y = Factors->value[G2::Y];

        double ** Z_X = Samples->value[G2::X];
        double ** Z_Y = Samples->value[G2::Y];

        double chol_y = sqrt(1.-rho*rho);

        size_t i;
        switch (nthreads){
            case 1:
                {
                    size_t path;
                    for (path = 0; path < npaths; ++path){
                        if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                        else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                        if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                        else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                    }
                    break;
                }
            default:
                {
                    size_t last_thread = nthreads-1, paths_per_thread = (npaths - (npaths % nthreads))/nthreads;

                    thread * threads = new thread[nthreads];
                    for (i = 0; i < last_thread; ++i){
                        threads[i] = thread([this,i,paths_per_thread,a,b,vol1,vol2,rho,chol_y,t,Z_X,Z_Y,X,Y]()->void{
                            size_t begin = i*paths_per_thread;
                            size_t end   = (i+1)*paths_per_thread;
                            size_t path;
                            for (path = begin; path < end; ++path){
                                if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                                else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                                if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                                else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                            }
                        });
                    }

                    threads[last_thread] = thread([this,npaths,last_thread,paths_per_thread,a,b,vol1,vol2,rho,chol_y,t,Z_X,Z_Y,X,Y]()->void{
                        size_t begin = last_thread*paths_per_thread;
                        size_t path;
                        for (path = begin; path < npaths; ++path){
                            if (a==0){ X[path] = vol1*sqrt(t)*Z_X[path][0]; }
                            else{      X[path] = vol1*sqrt(.5*(1.-exp(-2.*a*t))/a)*Z_X[path][0]; }
                            if (b==0){ Y[path] = vol2*sqrt(t)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                            else{      Y[path] = vol2*sqrt(.5*(1.-exp(-2.*b*t))/b)*(rho*Z_X[path][0]-chol_y*Z_Y[path][0]); }
                        }
                    });

                    for (i = 0; i < nthreads; ++i){
                        threads[i].join();
                    }

                    thread garbage_collection = thread([threads]()->void{ delete [] threads; }); 
                    garbage_collection.detach();
                    break;
                }// default
        }// switch(threads)
        clearDirty(G2::EVOLUTION);
    }//if (isDirty(G2::EVOLUTION))
    return Factors;
}

double G2PP::M(double x, double y, double T_t){
    double a = Coefs[G2::A];
    double b = Coefs[G2::B];
    
    double fact1 = x * (a ? (1.-exp(-a*T_t))/a : T_t);
    double fact2 = y * (b ? (1.-exp(-b*T_t))/b : T_t);

    return (fact1 + fact2);
}

double G2PP::V(double T_t){
    double a = Coefs[G2::A];
    double b = Coefs[G2::B];
    double vol1 = Coefs[G2::SIGMA_1];
    double vol2 = Coefs[G2::SIGMA_2];
    double rho = Coefs[G2::RHO];

    double fact1 = vol1*vol1*(a ? (T_t+2./a*exp(-a*T_t)-.5/a*exp(-2.*a*T_t)-1.5/a)/a/a : T_t);
    double fact2 = vol2*vol2*(b ? (T_t+2./a*exp(-b*T_t)-.5/b*exp(-2.*b*T_t)-1.5/b)/b/b : T_t);
    double fact12= 0;

    if ( a && b ){
        fact12 = vol1*vol2*(T_t + (exp(-a*T_t)-1.)/a + (exp(-b*T_t)-1.)/b - (exp(-(a+b)*T_t)-1.)/(a+b))/a/b;
    }else{
        if ( !a && !b ){
            fact12 = vol1*vol2*T_t;
        }else if (a && !b){
            fact12 = vol1*vol2*(T_t + (exp(-a*T_t)-1.)/a)/a;
        }else{
            fact12 = vol1*vol2*(T_t + (exp(-b*T_t)-1.)/b)/b;
        }
    }
    return (fact1 + fact2 + 2.*rho*fact12);
}
