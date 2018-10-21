#include "Optimization.h"

Optimization::Optimization(){}

Optimization::~Optimization(){
    if (Generator) {
        delete Generator;
    }
}


void calibrate(RateModel* model, RateInstrument* instruments, size_t num_instrs, double* weights, size_t max_iter, double precision, double k, size_t loss_trials){

}

void measureGradient(double * gradient, size_t num_params, RateModel * model, RateInstrument * instruments, size_t num_instrs, double * weights){

}

double loss_function (RateInstrument * instruments, size_t num_instrs, double * weights, int order = 2){
    double loss = 0;
    size_t i;
    for ( i = 0; i < num_instrs; ++i){
        loss += weights[i] * instruments[i].getLoss(order);
    }
    return loss;
}

double avg_loss (G2PP * model, RateInstrument * instruments, size_t num_instrs, double * weights, int ntrials, int order = 2){
    size_t i;
    double avg = 0;
    for ( i = 0; i < num_instrs; ++i){
        model->markDirtyAll(); // ensure simulation process kicks in to recalculate
        avg += fabs(loss_function(instruments, num_instrs, weights, order));
    }
    return (avg/(double)ntrials);
}
