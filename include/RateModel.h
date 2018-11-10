/*
 * ==========================================================================
 *
 *       Filename:  RateModel.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-24
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  James Ding
 *                  james.ding.illinois@gmail.com
 *
 * ==========================================================================
 */

#ifndef RATEMODEL_H
#define RATEMODEL_H

#include <string>
#include "Enums.h"

using namespace std;

struct RateModelError{
    int Code;
    string Message;

    size_t Line;
    string File;
    string Function;

    RateModelError(int code, string message): Code(code), Message(message){}
};

class RateModel{
    protected:
        RateModelType Type;
        string        Description;
    public:
        RateModelError ModelError;
        RateModel(RateModelType type, string description) : Type(type), Description(description), ModelError(1,"Success!") {}
        
        virtual ~RateModel(){

        }

        virtual RateModelType getModelType(){
            return Type;
        }

        virtual string getModelDescription(){
            return Description;
        }

        //Usage: getZCBP(double t, double T)
        //Output: P(t,T)
        virtual double getZCBP(double, double){
            return 0.0;
        }

        //Usage: getZCBP(double * prices, double * t, double T, size_t len)
        //Output: prices = [P(t[0],T) ,..., P(t[nterms-1], T)]
        virtual void getZCBP(double *, double *, double, size_t){ }

        //Usage: getZCBP(double * prices, double t, double * T, size_t len)
        //Output: prices = [P(t,T[0]) ,..., P(t, T[nterms-1])]
        virtual void getZCBP(double *, double, double *, size_t){ }

};

#endif
