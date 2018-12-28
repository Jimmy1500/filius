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
        WARN("precision must be greater than 0")
#endif
    }
    if ( k <= 0.0 ) {
        k = 0.01;
#ifdef __DEBUG__
        WARN("k must be greater than 0")
#endif
    }
    if ( alpha <= 0.0 ) {
        alpha = 0.01;
#ifdef __DEBUG__
        WARN("alpha must be greater than 0")
#endif
    }

    switch (model->getModelType()){
        case RMT_G2PP:
            if ( G2PP * g2pp = dynamic_cast<G2PP *>(model) ){
                const size_t num_params = 5;

                size_t keys[num_params] = {G2::A, G2::B, G2::SIGMA_1, G2::SIGMA_2, G2::RHO};
                double curr_guess[num_params] = {0., 0., 0., 0., 0.};
                double next_guess[num_params] = {0., 0., 0., 0., 0.};
                double best_guess[num_params] = {0., 0., 0., 0., 0.};
                double gradient[num_params]   = {0., 0., 0., 0., 0.};

                double factor = precision+0.1;
                double curr_temp = avg_loss(g2pp, instrs, weights, num_instrs, num_trials);
                double next_temp = curr_temp, best_temp = curr_temp, prob = 0.;

                g2pp->getParameters(keys, curr_guess, num_params);
                Generator = new Rand<double>((num_params+1),1,0,1);

                size_t i, iter = 0, optimum_reached = 0, tiny_improvement = 0;
                do{
                    if ( curr_temp <= precision ) { break; }

                    getGradient(gradient, keys, num_params, model, instrs, weights, num_instrs);
                    for (i = 0; i < num_params; ++i){
                        // stochastic descent (simulated annealing w. local optimizer)
                        next_guess[i] = curr_guess[i] - alpha * gradient[i] + Generator->nrand(0, i);
                    }
                    applyBoundaries(keys, next_guess, num_params);
                    g2pp->setParameters(keys, next_guess, num_params);
                    next_temp = avg_loss(g2pp, instrs, weights, num_instrs, num_trials);
#ifdef __DEBUG__
                    cout<<"### Iteration "<<(iter+1)<<" ###"<<endl;
                    cout<<"### Current gradient: ";
                    for (i = 0; i < num_params; ++i){
                        cout<<gradient[i]<<"; ";
                    }
                    cout<<endl;

                    cout<<"### Current configuration: ";
                    for ( i = 0; i < num_params; ++i ){
                        cout<<curr_guess[i]<<"; ";
                    }
                    cout<<endl;
                    cout<<"### Current temperature: "<<curr_temp<<endl;

                    cout<<"### Next configuration: ";
                    for ( i = 0; i < num_params; ++i ){
                        cout<<next_guess[i]<<"; ";
                    }
                    cout<<endl;
                    cout<<"### Next temperature: "<<next_temp<<endl;
#endif
                    if ( isnan(next_temp) ){
#ifdef __DEBUG__
                        DEBUG("### Adjacent temperature not measurable, trying again...")
#endif
                        continue;
                    } else if ( next_temp <= precision ) {
#ifdef __DEBUG__
                        DEBUG("### Adjacent temperature low enough, configuration is accepted as optimal solution!")
#endif
                        g2pp->getParameters(keys, best_guess, num_params);
                        curr_temp = next_temp; best_temp = next_temp;
                        break; //accept current state
                    } else if ( next_temp < curr_temp ){
                        g2pp->getParameters(keys, curr_guess, num_params);
                        curr_temp = next_temp;

                        if ( next_temp < best_temp ){
                            factor = 1 - next_temp/best_temp; //temperature improvement measure
                            g2pp->getParameters(keys, best_guess, num_params);
                            best_temp = next_temp;

                            if ( factor < precision ) {
                                tiny_improvement++;
                            }
                        }
#ifdef __DEBUG__
                        DEBUG("### Adjacent temperature is cooler (better), configuration accepted with improvment factor: "<<factor)
#endif
                    } else {
#ifdef __DEBUG__
                        DEBUG("### Adjacent temperature is hotter (worse), determining transition probability...")
#endif
                        prob = exp( (curr_temp - next_temp) / (k * curr_temp) );
#ifdef __DEBUG__
                        DEBUG("transition probability:  "<<prob)
#endif
                        if ( Generator->urand(0, num_params) <= prob ){
                            g2pp->getParameters(keys, curr_guess, num_params);
                            curr_temp = next_temp;
#ifdef __DEBUG__
                            DEBUG("transition accepted...")
#endif
                        } else{
                            g2pp->setParameters(keys, curr_guess, num_params);
#ifdef __DEBUG__
                            DEBUG("transition rejected...")
#endif
                        }
                    }

                    if ( isZero(gradient, num_params, precision) ){
                        ++optimum_reached;
                    } else {
                        optimum_reached = 0;
                    }
#ifdef __DEBUG__
                    DEBUG("local optimum reached: "<<optimum_reached)
                    cout<<endl;
#endif
                    ++iter;
                } while ( optimum_reached < 49 && tiny_improvement < 7 && iter < max_iter );

                if ( best_temp < curr_temp ){
                    g2pp->setParameters(keys, best_guess, num_params);
                }
#ifdef __DEBUG__
                DEBUG("### Accepted state ###");
                DEBUG("### Accepted temperature: "<<best_temp)
                DEBUG("### Accepted configuration: ")
                for ( i = 0; i < num_params-1; ++i ){
                    cout<<best_guess[i]<<",";
                }
                DEBUG(best_guess[num_params-1])
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
    size_t i;
    for ( i = 0; i < num_params; ++i ){
        if ( fabs(gradient[i]) > precision ){
            return false;
        }
    }
    return true;
}
