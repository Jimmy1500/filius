#include "Optimization.h"

Optimization::Optimization() : Generator(nullptr){}

Optimization::~Optimization(){
    if (Generator) {
        delete Generator;
    }
}


void Optimization::calibrate (RateModel* model, RateInstrument* instrs, double* weights, size_t num_instrs, size_t max_iter, double precision, double k, double alpha, size_t num_trials){
    if ( precision <= 0.0 ) {
        precision = 1.e-12;
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
                double gradient[num_params]   = {0., 0., 0., 0., 0.};

                double factor = precision+0.1;
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

                    for (i = 0; i < num_params; ++i){
                        // stochastic descent (simulated annealing w. local optimizer)
                        next_guess[i] = curr_guess[i] - alpha * gradient[i] + Generator->nrand(0, i);
                    }
                    applyBoundaries(keys, next_guess, num_params);
                    g2pp->setParameters(keys, next_guess, num_params);
#ifdef __DEBUG__
                    cout<<"### Iteration "<<(iter+1)<<" ###"<<endl;
                    cout<<"### Current gradient: ";
                    for (i = 0; i < num_params; ++i){
                        cout<<gradient[i]<<"; ";
                    }
                    cout<<endl;

                    cout<<"### Current temperature: "<<curr_temp<<endl;
                    cout<<"### Current configuration: ";
                    for ( i = 0; i < num_params; ++i ){
                        cout<<curr_guess[i]<<"; ";
                    }
                    cout<<endl;

                    next_temp = avg_loss(g2pp, instrs, weights, num_instrs, num_trials);
                    cout<<"### Next temperature: "<<next_temp<<endl;
                    cout<<"### Next configuration: ";
                    for ( i = 0; i < num_params; ++i ){
                        cout<<next_guess[i]<<"; ";
                    }
                    cout<<endl;
#endif
                    if ( isnan(next_temp) ){
                        cout<<"Adjacent temp not valid, trying again..."<<endl;
                        continue;
                    } else if ( next_temp <= precision ) {
                        cout<<"Adjacent temp low enough, configuration is optimal!"<<endl;
                        g2pp->getParameters(keys, curr_guess, num_params);
                        curr_temp = next_temp;
                        break; //accept current state
                    } else if ( next_temp < curr_temp ){
                        cout<<"Adjacent temp is better, configuration accepted..."<<endl;
                        g2pp->getParameters(keys, curr_guess, num_params);
                        factor = next_temp/initial_temp; //temperature improvement measure
                        curr_temp = next_temp;
                    } else {
                        cout<<"Adjacent temp is worse, determining transition probability..."<<endl;
                        prob = exp( (curr_temp - next_temp) / (k * curr_temp) );
                        cout<<"transition probability:  "<<prob<<", transition factor:  "<<factor<<", precision:  "<<precision<<endl;
                        if ( Generator->urand(0, num_params) <= prob ){
                            g2pp->getParameters(keys, curr_guess, num_params);
                            curr_temp = next_temp;
                            cout<<"transition accepted..."<<endl;
                        } else{
                            g2pp->setParameters(keys, curr_guess, num_params);
                            cout<<"transition rejected..."<<endl;
                        }
                    }
                    cout<<endl;
                    ++iter;
                }while ( !isZero(gradient, num_params, precision) && iter < max_iter && factor > precision );
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
    double delta = 1.e-9, f_left, f_right;

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
#ifdef __REGEN__
    size_t i;
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

bool Optimization::isZero(double * gradient, size_t num_params, double precision){
    double sum = 0.;
    size_t i;
    for ( i = 0; i < num_params; ++i ){
        sum += fabs(gradient[i]);
        if ( sum > precision ){
            return false;
        }
    }
    return true;
}
