/*
 * ==========================================================================
 *
 *       Filename:  Simulation.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-16
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
 *
 * ==========================================================================
 */


/* 
                      ******************
                     *                **
             ndims *                *  *
                 *     npaths     *    *
                ******************     *
                *                *     *
          nterms*                *    *
                *                *  *
                *                **
                ******************
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Matrix.h"
#include "Rand.h"
#include <thread>
#include <cmath>
//#include<iostream>

using namespace std;

class Simulation{
    private:
        double mu;
        double sigma;
        Rand<double> * Generator;

        void fitGenerator(size_t, size_t);
    public:
        Simulation();
        Simulation(double, double);
        ~Simulation();
        
        void setRandAttributes(double=0.0, double=1.0);
        /* match_dims = false(default): Dimensions of Random Generator (paths, dims) is [ 1, nthreads ]
         * match_dims = true: Dimensions of Random Generator (paths, dims) is the same as those of target random matrix  */
        void ugen(mat3d &, size_t, bool=false);//populate 3D matrix(double) with uniformly distributed number, i.i.d
        void ngen(mat3d &, size_t, bool=false);//populate 3D matrix(double) with uniformly distributed number, i.i.d
        
        Rand<double> * getGenerator();
};

#endif //SIMULATION_H
