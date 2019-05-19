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
        
        void Calibrate (RateModel *, RateInstrument *, double *, size_t, size_t, double precision=1.e-12, double k=0.01, double alpha=0.01, size_t num_trials=1);
        void GetGradient (double *, size_t *, size_t, RateModel *, RateInstrument *, double *, size_t) const;
        void ApplyBoundaries(size_t *, double *, size_t) const;
        bool IsZero(double *, size_t, double) const;

        double GetLoss (RateInstrument *, double *, size_t,  size_t order = 2) const;
        double GetAvgLoss (RateModel *, RateInstrument *, double *, size_t,  size_t , size_t order = 2) const;
};
#endif
