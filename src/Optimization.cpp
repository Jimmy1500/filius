#include "Optimization.h"

Optimization::Optimization(){}

Optimization::~Optimization(){
    if (Generator) {
        delete Generator;
    }
}


void Optimization::calibrate (RateModel* model, RateInstrument* instruments, size_t num_instrs, double* weights, size_t max_iter, double precision, double k, size_t loss_trials){

}

void Optimization::getGradient (double * gradient, size_t num_params, RateModel * model, RateInstrument * instruments, size_t num_instrs, double * weights){
    double f_left, f_right;
    double delta = 0.0001;

    size_t i;
    switch (model->getModelType()){
        case RMT_G2PP:
            if ( G2PP * g2pp = dynamic_cast<G2PP *>(model) ){
                double param_i;
                for ( i = 0; i < num_params; ++i){
                    param_i = g2pp->getParameter(i);
                    g2pp->setParameter(i, (param_i-delta));
                    f_left = loss_function(instruments, num_instrs, weights);

                    g2pp->setParameter(i, (param_i+delta));
                    f_right = loss_function(instruments, num_instrs, weights);

                    gradient[i] = (f_right - f_left) / delta*0.5;
                    g2pp->setParameter(i, param_i);
                }
            }
            break;
        case RMT_BLACK:
            break;
        default:
            break;
    }
}

double Optimization::loss_function (RateInstrument * instruments, size_t num_instrs, double * weights, size_t order){
    double loss = 0;
    size_t i;
    for ( i = 0; i < num_instrs; ++i){
        loss += weights[i] * instruments[i].getLoss(order);
    }
    return loss;
}

double Optimization::avg_loss (RateModel * model, RateInstrument * instruments, size_t num_instrs, double * weights, size_t num_trials, size_t order){
    size_t i;
    double avg = 0;
    for ( i = 0; i < num_trials; ++i){
        model->markDirtyAll(); // ensure simulation process kicks in to recalculate
        avg += fabs(loss_function(instruments, num_instrs, weights, order));
    }
    return (avg/(double)num_trials);
}
