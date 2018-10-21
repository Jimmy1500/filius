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
