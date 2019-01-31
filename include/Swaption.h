/*
 * ==========================================================================
 *
 *       Filename:  Swaption.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018-02-09
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
 *
 * ==========================================================================
 */

/* 
Swaption Term Structures:
      First Cashflow: Terms[0]                                                   Maturity/Last Cashflow: Terms[NumTerms-1]
Settle- ||                                                                         ||
ment   \||/<========================== Terms(NumTerms) ==========================>\||/
|  tau  \/  tau[0]    tau[1]     ...       ...       ...     tau[n-3]      tau[n-2]\/
t       T[0]      T[1]      T[2]      T[3]      T[4]      ...        T[n-2]        T[n-1]
*       *         *         *         *         *       *            *             *
--------------------------------------------------------------------------------------
Simulated ZCB Prices:
  | <============ P.npaths ============> |
  P(t,T[0])[0]    ...  P(t,T[0])[npaths-1]   ----
  P(t,T[1])[0]    ...  P(t,T[1])[npaths-1]    /\
       *           *         *               /||\
       *           *         *                ||
       *           *         *             P.nterms
       *           *         *                ||
       *           *         *               \||/
       *           *         *                \/
  P(t,T[n-1])[0]  ...  P(t,T[n-1])[npaths-1] ----
*/

#ifndef SWAPTION_H
#define SWAPTION_H

#include "Matrix.h"
#include "G2PP.h"
#include "RateInstrument.h"

#ifndef MAX
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#endif

class Swaption : public RateInstrument{
    private:
        RateModel* Model;

        double* Terms;
        size_t  NumTerms;

        size_t *Params;       //Parameters: NTL, STK, STL
        size_t *Flags;        //Dirty marker for Params
        volatile size_t Dirty;

        mat2d  *Prices;       //Intermediate Simulated bond price
        double  Payoff;       //Intermediate Payoff

    public:
        Swaption();
        Swaption(RateModel *);
        Swaption(Swaption &);
        ~Swaption();

        double getModelValue();
        double getModelValue(double, double, double, double *, size_t);

        //--------------------------------Setters----------------------------------------------
        inline void setRateModel(RateModel * model){
            if (model){
                if (Model){
                    Description.erase(14);
                    if (Model->getModelDescription() != model->getModelDescription()){
                        Description.append(model->getModelDescription());
                    }
                }else{ Description.append(model->getModelDescription()); }
                Model = model;
            }
        }

        inline void setParameter(size_t key, double value){
            Params[key] = value;
            markDirtyFrom(Flags[key]);
        }

        inline void setParameters(size_t * keys, double * values, size_t len){
            size_t i;
            if (len > SWPT::NUM_PARAMS){
                REPORT_ERROR(Model->ModelError, "Too many parameters to be expected",0)
            }
            for (i=0; i<len; ++i){ setParameter(keys[i], values[i]); }
        }

        inline void setTerms(double * terms, size_t nterms) { Terms = terms; NumTerms = nterms; markDirtyAll(); }

        //--------------------------------Getters-----------------------------------------------
        constexpr RateModel* getModel()            const    { return Model; }
        constexpr double* getTerms()               const    { return Terms; }
        constexpr size_t  getNumTerms()            const    { return NumTerms; }
        constexpr double  getParameter(size_t key) const    { return Params[key]; }

        //--------------------------------Utilities-----------------------------------------------
        constexpr double maxDiff(double float_leg, double fix_leg) const {
            return ( (float_leg > fix_leg) ? (float_leg - fix_leg) : (0.) );
        }

        //-----------Procedural Step Components Recalculation Bitwise Marker/Verifier------------
        inline void markDirtyFrom(size_t step){
            if (step){
                size_t i;
                for (i=step; i<SWPT::NUM_PROCS; i<<=1u){
                    Dirty|=i;
                }
            }
        }

        inline void markDirty(size_t step){
            Dirty|=step;
        }

        inline void markDirtyAll(){
            size_t i;
            for (i=SWPT::GET_ZCBP; i<SWPT::NUM_PROCS; i<<=1u){
                Dirty|=i;
            }
            if (Model){
                switch (Model->getModelType()){
                    case RMT_G2PP:
                        if ( G2PP * g2pp = dynamic_cast<G2PP *>(Model) ){
                            g2pp->markDirtyAll();
                        }else{
                            REPORT_ERROR(Model->ModelError, "Unable to resolve rate model type (G2PP)", 0)
                        }
                        break;
                    case RMT_BLACK:
                        break;
                    default:
                        break;
                }
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

        inline size_t isDirty(){
            return ( Dirty > 0 );
        }
};

#endif
