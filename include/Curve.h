/*
 * ==========================================================================
 *
 *       Filename:  Curve.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-17
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
 *
 * ==========================================================================
 */
#ifndef CURVE_H
#define CURVE_H

#include <cstddef>
#include <cmath>

using namespace std;

#define DeepCopy_YC(terms, prices)      \
        size_t i;                       \
        for (i = 0; i < Length; ++i){   \
            Terms[i] = terms[i];        \
            Values[i] = prices[i];      \
        }

#define ShallowCopy_YC(terms, prices)   \
        Terms = terms;                  \
        Values = prices;

class Curve{
    private:
        double *Terms;     //Terms
        double *Values;    //Prices w.r.t settlement(Time), i.e. { P(Time, Terms[i]) }
        size_t  Length;    //Length of the yield curve
    public:
        Curve(Curve *);
        Curve(double*, double*, size_t);
        Curve(double, double*, double*, size_t);
        ~Curve();

        //Utilities:
        //ORDED BY Terms from least to most
        void interpolateTerms(double, size_t&, int&); //Outputs the nearest 2 indices for inter/extrapolation
        double P(double);                             //Log-Linear Interpolation
        double P(double, double);                     //P(t,T) = P(0,T)/P(0,t)

        //Getters/Setters
        constexpr double * getTerms()  const;
        constexpr double * getValues() const;
        constexpr size_t   getLength() const;

};

#endif
