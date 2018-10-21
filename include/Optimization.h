/*
 * ==========================================================================
 *
 *       Filename:  Optimization.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018-10-21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jimmy1500
 *                  lighteningmagic@gmail.com
 *
 * ==========================================================================
 */
#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H
#include <math.h>
#include <iostream>
#include "Rand.h"
#include "G2PP.h"
#include "Swaption.h"

class Optimization{
    private:
        Rand<double> * Generator;
    public:
        Optimization();
        ~Optimization();
        
        void calibrate (RateModel *, RateInstrument *, size_t, double *, size_t, double, double k=0.01, size_t loss_trials=1);
        void getGradient (double *, size_t, RateModel *, RateInstrument *, size_t, double *);
        double loss_function (RateInstrument * , size_t , double * , size_t order = 2);
        double avg_loss (RateModel * , RateInstrument * , size_t , double * , size_t , size_t order = 2);
};
#endif
