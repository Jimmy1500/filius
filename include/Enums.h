/*
 * ==========================================================================
 *
 *       Filename:  G2Enums.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-26
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
 *
 * ==========================================================================
 */
#ifndef ENUMS_H
#define ENUMS_H

#include <mutex>

static mutex mtx;

enum RateModelType {RMT_G2PP=0, RMT_BLACK, NUM_MODEL_TYPES};

namespace G2
{
    enum PROC{
        GENERATION=1,
        EVOLUTION=2,
        SIMULATION=4, 
        CALCULATION=8,
        NUM_PROCS
    };

    enum PARAM{
        A=0, 
        B, 
        SIGMA_1, 
        SIGMA_2, 
        RHO,
        PC_A,      //predictor corrector alpha (only used when nterms > 1)
        PC_B,      //predictor corrector eta   (only used when nterms > 1)
        NUM_PARAMS
    };

    //Peripheries
    enum PERI{
        NTERMS=0, //nterms = 1: analytical solution(s) exist for underlying model SDE(s); nterms > 1: predictor corrector/euler maruyama
        NPATHS,
        NDIMS,
        NTHREADS,
        NUM_PERIS
    };

    enum DIM{
        X=0, 
        Y, 
        NUM_DIMS
    };
}

namespace SWPT
{
    enum PROC{
        GET_ZCBP=1,
        GET_PAYOFF=2,
        GET_VALUE=4,
        NUM_PROCS
    };

    enum PARAM{
        NTL=0,
        STK,
        STL,
        NUM_PARAMS
    };
}

#endif
