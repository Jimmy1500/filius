#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>

#include "Optimization.h"

using namespace std;

int main(int argc, char *argv[]) {
    RateModel * model = new G2PP();
    cout<<model->getModelDescription()<<endl;
    G2PP * g2pp = static_cast<G2PP*>(model);

    double t = 2; double T = 30;
    const int N = 11;
    double stl[N]={1./12., 0.25, .5, 1., 2., 3., 5., 7., 10., 20., 30.};
    double price [N]={
        exp(-.0193*stl[0]), exp(-.0203*stl[1]), exp(-.0222*stl[2]), 
        exp(-.0245*stl[3]), exp(-.0267*stl[4]), exp(-.0278*stl[5]), 
        exp(-.0287*stl[6]), exp(-.0296*stl[7]), exp(-.0300*stl[8]),
        exp(-.0307*stl[9]), exp(-.0313*stl[10])
    };
    Curve YieldCurve = Curve(stl, price, N);
    RateInstrument * instr = new Swaption(model);
    Swaption * swaption = dynamic_cast< Swaption * > (instr);
    
    try{
        size_t peris_keys[] = {G2::NTERMS, G2::NPATHS, G2::NDIMS, G2::NTHREADS};
        size_t peris[] = {1,777777,2,16};
        size_t coefs_keys[] = {G2::A, G2::B, G2::SIGMA_1, G2::SIGMA_2, G2::RHO, G2::PC_A, G2::PC_B};
        double coefs[] = {0.10, 0.10, 0.042, 0.045, -0.89, 0.8, 0.5};
        
        g2pp->setYieldCurve(&YieldCurve);
        g2pp->setPeripheries(peris_keys, peris, G2::PERI::NUM_PERIS);
        g2pp->setParameters(coefs_keys, coefs, G2::PARAM::NUM_PARAMS);
        
        size_t i,j,nterms = 7;
        double Time, AVG;
        double PVs[nterms], Stls[nterms], Mats[nterms];
        cout.precision(10);
        cout<<"Current:P("<<t<<","<<T<<"): "<<YieldCurve.P(t,T)<<endl;

        //getZCBP(t,T) -> P(t,T)
        cout<<"Testing getZCBP(t,T)"<<endl;
        Time = 0.0; AVG = 0.0;
        for (i = 0; i < nterms; i++){
            auto start = chrono::high_resolution_clock::now();
            double PV = g2pp->getZCBP(t,T);
            auto stop = chrono::high_resolution_clock::now();
            chrono::duration<double, milli> elapsed = stop-start;

            cout<<"Future::P("<<t<<","<<T<<"): "<<PV<<endl;
            Time += elapsed.count(); AVG += PV;
            g2pp->clearSimulation();
        }
        cout<<"# of trials:\t"<<i<<endl;
        cout<<"Average cost:\t"<<Time/(double)i<<" milliseconds"<<endl;
        cout<<"Average value:\t"<<AVG/(double)i<<endl<<endl;

        //getZCBP(t[nterms],T) -> P(t[0],T), ..., P(t[nterms-1],T)
        cout<<"Testing getZCBP(double *prices,double *t,T,nterms)"<<endl;
        Time = 0.0; AVG = 0.0;
        fill(PVs, PVs+nterms, 0.0); fill(Stls, Stls+nterms, t);
        //for (i = 1; i < nterms; i++){ Stls[i] += i*.01; }
        auto start = chrono::high_resolution_clock::now();
        g2pp->getZCBP(PVs, Stls, T, nterms);
        auto stop = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> elapsed = stop-start;

        Time = elapsed.count();
        for (i = 0; i < nterms; i++){
            AVG += PVs[i];
            cout<<"Future::P("<<Stls[i]<<","<<T<<"): "<<PVs[i]<<endl;
        }
        cout<<"# of trials:\t"<<i<<endl;
        cout<<"Average cost:\t"<<Time/(double)i<<" milliseconds"<<endl;
        cout<<"Average value:\t"<<AVG/(double)i<<endl<<endl;

        g2pp->clearSimulation();
        //getZCBP(t,T[nterms]) -> P(t,T[0]), ..., P(t,T[nterms-1])
        cout<<"Testing getZCBP(double *prices,t,double *T,nterms)"<<endl;
        Time = 0.0; AVG = 0.0;
        fill(PVs, PVs+nterms, 0.0); fill(Mats, Mats+nterms, T);
        start = chrono::high_resolution_clock::now();
        g2pp->getZCBP(PVs, t, Mats, nterms);
        stop = chrono::high_resolution_clock::now();
        elapsed = stop-start;

        Time = elapsed.count();
        for (i = 0; i < nterms; i++){
            AVG += PVs[i];
            cout<<"Future::P("<<t<<","<<Mats[i]<<"): "<<PVs[i]<<endl;
        }
        cout<<"# of trials:\t"<<i<<endl;
        cout<<"Average cost:\t"<<Time/(double)i<<" milliseconds"<<endl;
        cout<<"Average value:\t"<<AVG/(double)i<<endl<<endl;
        g2pp->clearSimulation();

        cout<<"Testing getZCBP(double **prices,t,double *T,nterms,npaths)"<<endl;
        AVG = 0.0;
        mat2d prices(nterms, g2pp->getPeriphery(G2::PERI::NPATHS));
        start = chrono::high_resolution_clock::now();
        g2pp->getZCBP(prices.value, t, Mats, prices.nterms, prices.npaths);
        stop = chrono::high_resolution_clock::now();
        elapsed = stop-start;
        Time = elapsed.count();

        for (j = 0; j < prices.nterms; j++){
            PVs[j] = 0.0;
            for (i = 0 ; i < prices.npaths; i++){
                PVs[j] += prices.value[i][j];
            }
            PVs[j] /= (double)prices.npaths;
            AVG += PVs[j];
            cout<<"Future::P("<<t<<","<<Mats[j]<<"): "<<PVs[j]<<endl;
        }
        cout<<"# of trials:\t"<<j<<endl;
        cout<<"Average cost:\t"<<Time/(double)j<<" milliseconds"<<endl;
        cout<<"Average value:\t"<<AVG/((double)j)<<endl<<endl;

        const size_t NN = 12, num_params = 3;
        double terms[NN] = {t+3./12, t+6./12, t+9./12, t+12./12, t+15./12., t+18./12, t+21./12, t+24./12, t+27./12, t+30./12, t+33./12, t+36./12};
        double notional = 100., strike = .02;
        size_t swpt_keys  [num_params] = {SWPT::NTL, SWPT::STK, SWPT::STL};
        double swpt_values[num_params] = { notional,    strike,       t  };

        swaption->setParameters(swpt_keys, swpt_values, num_params);
        swaption->setTerms(terms, NN);
        swaption->setMarketValue(8.01);

        cout<<"Instrument Type: "<<(swaption->getInstrumentType()==RIT_SWAPTION ? "Swaption" : "Unknown")<<endl;
        cout<<"Pricing "<<swaption->getInstrumentDescription()<<endl;

        //g2pp->setPeriphery(G2::PERI::NTHREADS, 1);
        Time = 0.0; AVG = 0.0;
        for (i = 0; i < nterms; i++){
            auto start = chrono::high_resolution_clock::now();
            double swaption_model_value = swaption->getModelValue();
            cout<<"Swaption(NTL:"<<notional<<",STK:"<<strike<<",STL:"<<t<<",MAT:"<<terms[NN-1]<<",FRQ:"<<terms[1]-terms[0]<<",NTERMS:"<<NN<<"): "<<swaption_model_value<<endl;
            auto stop = chrono::high_resolution_clock::now();
            chrono::duration<double, milli> elapsed = stop-start;
            Time += elapsed.count(); AVG += swaption_model_value;
        }
        cout<<"Average Cost: "<<Time/(double)i<<" milliseconds"<<endl;
        cout<<"Average Price: "<<AVG/(double)i<<endl<<endl;;

        const size_t num_instrs = 1, max_iter = 212;
        Swaption swpts[num_instrs] = {*swaption};
        double weights[num_instrs] = {1.0};

        for (i = 0; i < num_instrs; ++i){
            cout<<"Setting swaption with allocated weight: "<<weights[i]<<endl;
            swpts[i].setParameters(swpt_keys, swpt_values, num_params);
            swpts[i].setTerms(terms, NN);
            swpts[i].setMarketValue(7.81 + .01*(double)i);
        }
        Optimization * opt = new Optimization();
        opt->calibrate(g2pp, swpts, weights, num_instrs, max_iter);
        delete opt;

        cout<<"### Target instrument model value: "<<endl;
        cout<<"### Swaption(NTL:"<<notional<<",STK:"<<strike<<",STL:"<<t<<",MAT:"<<terms[NN-1]<<",FRQ:"<<terms[1]-terms[0]<<",NTERMS:"<<NN<<"): "<<swpts[0].getMarketValue()<<endl;
        cout<<"### Accepted instrument model value: "<<endl;
        cout<<"### Swaption(NTL:"<<notional<<",STK:"<<strike<<",STL:"<<t<<",MAT:"<<terms[NN-1]<<",FRQ:"<<terms[1]-terms[0]<<",NTERMS:"<<NN<<"): "<<swpts[0].getModelValue()<<endl;

    }catch(int code){
        cout<< model->ModelError.Message<<"@"<<model->ModelError.Function<<":"<<model->ModelError.File<<":"<<model->ModelError.Line <<endl;
    }

    thread garbage = thread([instr, model]()->void{ delete instr; delete model; });
    garbage.detach();

    return 0;
}
