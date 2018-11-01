#include "Optimization.h"

Optimization::Optimization() : Generator(nullptr){}

Optimization::~Optimization(){
    if (Generator) {
        delete Generator;
    }
}


void Optimization::calibrate (RateModel* model, RateInstrument* instrs, double* weights, size_t num_instrs, size_t max_iter, double precision, double k, double alpha, size_t num_trials){
    if ( precision <= 0.0 ) {
        precision = 1.e-6;
#ifdef __DEBUG__
        DEBUG("precision must be greater than 0")
#endif
    }
    if ( k <= 0.0 ) {
        k = 0.01;
#ifdef __DEBUG__
        DEBUG("k must be greater than 0")
#endif
    }
    if ( alpha <= 0.0 ) {
        alpha = 0.01;
#ifdef __DEBUG__
        DEBUG("alpha must be greater than 0")
#endif
    }

    switch (model->getModelType()){
        case RMT_G2PP:
            if ( G2PP * g2pp = dynamic_cast<G2PP *>(model) ){
                const size_t num_params = 5;

                size_t keys[num_params] = {G2::A, G2::B, G2::SIGMA_1, G2::SIGMA_2, G2::RHO};
                double curr_guess[num_params] = {0., 0., 0., 0., 0.};
                double next_guess[num_params] = {0., 0., 0., 0., 0.};
                double gradient[num_params] = {0., 0., 0., 0., 0.};

                double factor = precision;
                double initial_temp = avg_loss(g2pp, instrs, weights, num_instrs, num_trials);
                double curr_temp = initial_temp, next_temp = initial_temp, prob = 0.;

                g2pp->getParameters(keys, curr_guess, num_params);
                Generator = new Rand<double>((num_params+1),1,0,1);

                size_t i, iter = 0;
                do{
                    if ( curr_temp <= precision ) {
                        break; //accept current state
                    } 

                    getGradient(gradient, keys, num_params, model, instrs, weights, num_instrs);
#ifdef __DEBUG__
                    cout<<"### iteration "<<(iter+1)<<" ###"<<endl;
                    cout<<"Current gradient: ";
                    for (i = 0; i < num_params; ++i){ 
                        cout<<gradient[i]<<",";
                    }
                    cout<<endl;
#endif
                    for (i = 0; i < num_params; ++i){ 
                        // stochastic descent (simulated annealing w. local optimizer)
                        next_guess[i] += Generator->nrand(0, i) - alpha * gradient[i];
                    }
                    applyBoundaries(keys, next_guess, num_params);
                    g2pp->setParameters(keys, next_guess, num_params);

                    next_temp = avg_loss(g2pp, instrs, weights, num_instrs, num_trials);
                    if ( isnan(next_temp) ){
                        continue;
                    } else if ( next_temp <= precision ) {
                        curr_temp = next_temp;
                        g2pp->getParameters(keys, curr_guess, num_params);
                        break; //accept current state
                    } else if ( next_temp < curr_temp ){
                        curr_temp = next_temp;
                        g2pp->getParameters(keys, curr_guess, num_params);
                        factor = curr_temp/initial_temp; //temperature improvement measure
                    } else {
                        prob = exp( (curr_temp - next_temp) / (k * curr_temp) );
                        if ( Generator->urand(0, num_params) <= prob ){
                            curr_temp = next_temp;
                            g2pp->getParameters(keys, curr_guess, num_params);
                        } else{
                            g2pp->setParameters(keys, curr_guess, num_params);
                        }
                    }
                    ++iter;
#ifdef __DEBUG__
                    cout<<"Current temperature: "<<curr_temp<<endl;
                    cout<<"Next temperature: "<<next_temp<<endl;
                    cout<<"Current parametric configuration: ";
                    for ( i = 0; i < num_params-1; ++i ){
                        cout<<curr_guess[i]<<",";
                    }
                    cout<<curr_guess[num_params-1]<<endl;
                    cout<<"### Continue iteration ..."<<endl;
#endif
                }while ( iter < max_iter && factor >= precision );
#ifdef __DEBUG__
                cout<<"### Accepted state ###"<<endl;
                cout<<"Accepted temperature: "<<curr_temp<<endl;
                cout<<"Accepted parametric configuration: ";
                for ( i = 0; i < num_params-1; ++i ){
                    cout<<curr_guess[i]<<",";
                }
                cout<<curr_guess[num_params-1]<<endl;
#endif
            }
            break;
        case RMT_BLACK:
            break;
        default:
            break;
    }
}

void Optimization::getGradient (double * gradient, size_t * param_keys, size_t num_params, RateModel * model, RateInstrument * instrs, double * weights, size_t num_instrs){
    double f_left, f_right;
    double delta = 1.e-5;

    size_t i;
    switch (model->getModelType()){
        case RMT_G2PP:
            if ( G2PP * g2pp = dynamic_cast<G2PP *>(model) ){
                double param;
                for ( i = 0; i < num_params; ++i){
                    param = g2pp->getParameter(param_keys[i]);
                    g2pp->setParameter(param_keys[i], (param-delta));
                    f_left = loss_function(instrs, weights, num_instrs);

                    g2pp->setParameter(param_keys[i], (param+delta));
                    f_right = loss_function(instrs, weights, num_instrs);

                    gradient[i] = (f_right - f_left) / delta*0.5;
                    g2pp->setParameter(param_keys[i], param);
                }
            }
            break;
        case RMT_BLACK:
            break;
        default:
            break;
    }
}

double Optimization::loss_function (RateInstrument * instrs, double * weights, size_t num_instrs, size_t order){
    double loss = 0;
    size_t i;
    for ( i = 0; i < num_instrs; ++i ){
        loss += weights[i] * instrs[i].getLoss(order);
    }
    return loss;
}

double Optimization::avg_loss (RateModel * model, RateInstrument * instrs, double * weights, size_t num_instrs, size_t num_trials, size_t order){
    size_t i;
#ifdef __REGEN__
    double avg = 0.;
    for ( i = 0; i < num_trials; ++i ){
        model->markDirtyAll(); 
        avg += fabs(loss_function(instrs, weights, num_instrs, order));
    }
    return (avg/(double)num_trials);
#else
    return fabs(loss_function(instrs, weights, num_instrs, order));
#endif
}

void Optimization::applyBoundaries(size_t * keys, double * values, size_t num_params){
    size_t i;
    for ( i = 0; i < num_params; ++i) {
        switch (keys[i]){
            case G2::A:
            case G2::B:
            case G2::SIGMA_1:
            case G2::SIGMA_2:
                values[i] = ( values[i] < 0.0 ? 0.0 : values[i] );
                break;
            case G2::RHO:
                if ( values[i] > 1.0 ) {
                    values[i] = 1.0;
                }else if ( values[i] < 0.0 ) {
                    values[i] = 0.0;
                }
                break;
        }
    }
}
