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
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
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

#define DEEPCOPY_ARRAY(target, value, length)   \
        for (i = 0; i < length; ++i){           \
            target[i] = value[i];               \
        }

class Optimization{
    private:
        Rand<double> * Generator;
    public:
        Optimization();
        ~Optimization();
        
        void calibrate (RateModel *, RateInstrument *, size_t, double *, size_t, double, double k=0.01, double alpha=0.5, size_t loss_trials=3);
        void getGradient (double *, size_t *, size_t, RateModel *, RateInstrument *, double *, size_t);
        double loss_function (RateInstrument *, double *, size_t,  size_t order = 2);
        double avg_loss (RateModel *, RateInstrument *, double *, size_t,  size_t , size_t order = 2);
};
#endif
