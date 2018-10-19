/*
 * ==========================================================================
 *
 *       Filename:  G2++.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-25
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jimmy1500
 *                  lighteningmagic@gmail.com
 *
 * ==========================================================================
 */

/* 
G2++:                   r[t] = X[t] + Y[t] + Phi[t], X[0] = 0; Y[0] = 0;
                    where:                            
                        dX[t] = -a*X[t] + sigma1*dW[t] = -a*X[t] + sigma1*sqrt(dt)*z1[t]
                        dY[t] = -b*Y[t] + sigma2*dW[t] = -b*Y[t] + sigma2*sqrt(dt)*z2[t]
                        Phi[t]= f_m_(0, t) + sigma1^2/(2*a^2)*[1-e^(-a*t)]^2
                                           + sigma2^2/(2*b^2)*[1-e^(-b*t)]^2
                                           + rho*sigma1*sigma2/(a*b)*[1-e^(-a*t)]*[1-e^(-b*t)]
                        Integral(0,t)->Phi[u]du = Pm(0,t)*exp[-0.5*V(0,T)]
                        Integral(t,T)->Phi[u]du = Pm(0,T)/Pm(0,t)*exp{-0.5*[V(0,T)-V(0,t)]}

                        Integral(t,T)->(X[u] + Y[u])du ~ N [M(t,T), V(t,T)]
                    where:
                        M(t,T) = X[t]*{1-exp[-a*(T-t)]}/a + Y[t]*{1-exp[-b*(T-t)]}/b
                        V(t,T) = {T-t + exp[-a*(T-t)]*2/a - exp[-2*a*(T-t)]/(2*a)-3/(2*a)}*sigma1^2/a^2
                               + {T-t + exp[-b*(T-t)]*2/b - exp[-2*b*(T-t)]/(2*b)-3/(2*b)}*sigma2^2/b^2
                               + <T-t + {exp[-a*(T-t)]-1}/a - {exp[-b*(T-t)]-1}/b + {exp[-(a+b)*(T-t)]-1}/(a+b)>*2*rho*sigma1*sigma2/(a*b)

                    Theoretical Bond Price at settlement t w/ maturity T
                        P(t,T) = Pm(0,T)/Pm(0,t) * exp(0.5*[V(t,T)-V(0,T)+V(0,t)] - M(t,T))
-------------------------------------------------------------------------------------------
Simulation Data Structure Visualization:

Orstein-Uhlenbeck Speed:  [ a         ,      b     ]
                            ||               ||
                            ||   applys to   ||
                            ||               ||
                Matrix[dim=0][path=0][t]  Matrix[dim=1][path=0][t]
                            ||               ||
                           \||/    **       \||/    **
                            \/   *  *        \/   * <=======thread[i+1]
                           X[t]*<==========>Y[t]* *  *
                              *  *  *          * <======thread[i] 
Dynamics:                     *  *  *          *  *  *
                         term *  * *      term *  * *
                              *  *             *  *
                              **               **
                             /\               /\
                            /||\   **        /||\   **
                             ||  *  *         ||  * <=======thread[i+1]
                          z1[t]*<===rho===>z2[t]* *  *
                              *  *  *          * <======thread[i]
Samples(rho):                 *  *  *          *  *  *
                              *  * *           *  * *
                              *  *             *  *
                              **               **
                             /\               /\
                            /||\  applys to  /||\
                             ||               ||
Volatilies:               [ sigma1    ,      sigma2 ]


*/

#ifndef G2PP_H
#define G2PP_H

#include "Simulation.h"
#include "Curve.h"
#include "RateModel.h"
#include <mutex>

#ifdef __DEBUG__
extern mutex mtx;
#include <iostream>
#define DEBUG(MSG) mtx.lock(); cout<<"Warning: "<<MSG<<endl; mtx.unlock();
#endif

#ifndef REPORT_ERROR
#define REPORT_ERROR(Error, MSG,CODE)          \
        Error.Message = MSG;                   \
        Error.Code = CODE;                     \
        Error.Function = __func__;             \
        Error.Line = __LINE__;                 \
        Error.File = __FILE__;                 \
        throw Error.Code;
#endif

#ifndef M_XY
#define M_XY(X,Y,T_t)                          \
        (                                      \
            X*(a?(1.-exp(-a*T_t))/a:T_t)       \
                          +                    \
            Y*(b?(1.-exp(-b*T_t))/b:T_t)       \
        )
#endif

#ifndef V_XY
#define V_XY(T_t)                                                                           \
        (                                                                                   \
            vol1*vol1*(a ? (T_t+2./a*exp(-a*T_t)-.5/a*exp(-2.*a*T_t)-1.5/a)/a/a : T_t)      \
                                                +                                           \
            vol2*vol2*(b ? (T_t+2./a*exp(-b*T_t)-.5/b*exp(-2.*b*T_t)-1.5/b)/b/b : T_t)      \
                                                +                                           \
            2.*rho*vol1*vol2*                                                               \
            (                                                                               \
                a && b ? (T_t+(exp(-a*T_t)-1.)/a+(exp(-b*T_t)-1.)/b                         \
                                                -                                           \
                                (exp(-(a+b)*T_t)-1.)/(a+b))/a/b :                           \
                (                                                                           \
                    !a && !b ? T_t :                                                        \
                    (                                                                       \
                        a && !b ? (T_t+(exp(-a*T_t)-1.)/a)/a : (T_t+(exp(-b*T_t)-1.)/b)/b   \
                    )                                                                       \
                )                                                                           \
            )                                                                               \
        )
#endif

#ifndef SIM_ZCB
#define SIM_ZCB(YieldCurve,X,Y,t,T,T_t)                                                     \
        (                                                                                   \
                .5 * YieldCurve->P(t,T)                                                     \
                            *                                                               \
            exp( .5*( V_XY(T_t) - V_XY(T) + V_XY(t)) )                                      \
                            *                                                               \
          ( exp(-M_XY(X,Y,T_t)) + exp(-M_XY(-X,-Y,T_t)) )                                   \
        )
#endif

class G2PP : public RateModel{
    private:
        Simulation *Sim;        //Simulation Core, responsible for generating random numbers in bulk
        Curve      *YieldCurve; //Yield curve, contains current market condition, and inter/extrapolation utilities
        mat3d      *Samples;    //Random matrix containing: dimension 1: z1~N(0,1), i.i.d, dimension 2: z2~N(0,1), i.i.d
        mat2d      *Factors;    //Simulated factors X[t],Y[t], and X_[t],Y_[t] (antithetic)
        mat1d      *Results;    //Intermediate simulation results exp(-M(t,T))
        double      Settlement;
        double      Maturity;
        double      ZCB;        //ZCB Price P(t,T) previous calculated 

        volatile size_t Dirty;  //Binary process recalcuation flag
        size_t     *Flags;      //Procedure dirty flags associated with Model Coefficients, See G2Enum.h for more detail

    protected:
        //Model parameters 
        double *Coefs;          //Model Coefficients (a, b, vol1, vol2, etc.)
        size_t *Peris;          //Simulation Peripheries (ndims, npaths, nterms, nthreads)

    public:
        //Constructor: pass in terms, npaths, ndims + distributional paramters of simulation
        G2PP();
        G2PP(Curve *);
        ~G2PP();

        //----Calculation Services-----
        double getZCBP(double  , double);                                       //output expected ZCB price   P(t,T)
        void   getZCBP(double *, double *, double  , size_t);                   //output expected ZCB prices [P(t[0],T), ..., P(t[n-1],T)]
        void   getZCBP(double *, double  , double *, size_t);                   //output expected ZCB prices [P(t,T[0]), ..., P(t,T[n-1])]
        void   getZCBP(double **,double  , double *, size_t, size_t);           //output simulated ZCB prices
                                                                                // [ P(Xt_0,Yt_0,t,T[0]),        ..., P(Xt_npaths,Yt_npaths,t,T[0])         ]
                                                                                // [        ...                  ...,           ...                         ]
                                                                                // [ P(Xt_0,Yt_0,t,T[nterms-1]), ..., P(Xt_npaths,Xt_npaths,t,T[nterms-1])  ]

        mat2d* getFactors(double);                                              //output simulated factors
                                                                                // [    X(t)(0)    ,    Y(t)(0)     ]
                                                                                // [      ...             ...       ]
                                                                                // [ X(t)(npaths-1), Y(t)(npaths-1) ]

        double M(double, double, double);
        double V(double);

        //----Getters & Setters-----
        inline void setParameter(size_t key, double value){
            Coefs[key] = value;
            markDirtyFrom(Flags[key]);
        }

        inline void setParameters(size_t * keys, double * values, size_t len){
            size_t i;
            if (len > G2::NUM_PARAMS){
                REPORT_ERROR(ModelError, "Number of parameters too great to be possible",0)
            }
            for (i=0; i<len; ++i){
                setParameter(keys[i], values[i]);
            }
        }

        inline void setPeriphery(size_t key, size_t value){
            if (Peris[key] != value){
                Peris[key] = value;
                markDirtyFrom(G2::GENERATION);
            }
        }

        inline void setPeripheries(size_t * keys, size_t * values, size_t len){
            size_t i;
            for (i=0; i<len; ++i){
                setPeriphery(keys[i], values[i]);
            }
        }

        inline void setYieldCurve(Curve * curve){ //Shallow Copy
            this->YieldCurve = curve;
            markDirtyFrom(G2::CALCULATION);
        }

        inline double getParameter(size_t key){
            return Coefs[key];
        }

        inline size_t getPeriphery(size_t key){
            return Peris[key];
        }

        inline Simulation* getSimEngine(){
            return Sim;
        }

        //----Utility----
        inline void clearSimulation(){
            delete Samples;Samples = nullptr;
            delete Factors;Factors = nullptr;
            delete Results;Results = nullptr;
            markDirtyAll();
        }

        //-----Simulation Procedural Step Components Recalculation Bitwise Marker/Verifier-----
        inline void markDirtyFrom(size_t step){
            if (step){
                size_t i;
                for (i=step; i<G2::NUM_PROCS; i<<=1u){
                    Dirty|=i;
                }
            }
        }

        inline void markDirty(size_t step){
            Dirty|=step;
        }

        inline void markDirtyAll(){
            size_t i;
            for (i=G2::GENERATION; i<G2::NUM_PROCS; i<<=1u){
                Dirty|=i;
            }
        }

        inline void clearDirty(size_t step){
            Dirty&=~step;
        }

        inline void clearDirtyAll(){
            Dirty=0;
        }

        inline size_t isDirty(size_t step){
            return Dirty & step;
        }
};

#endif