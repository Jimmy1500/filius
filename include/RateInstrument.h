/*
 * ==========================================================================
 *
 *       Filename:  RateInstrument.h
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
#ifndef RATEINSTRUMENT_H
#define RATEINSTRUMENT_H

#include <string>

using namespace std;

enum RateInstrumentType{RIT_SWAPTION=0};

class RateInstrument{
    protected:
        RateInstrumentType Type;
        string Description;
        double MarketValue;
    public:
        RateInstrument(RateInstrumentType type,string description) : Type(type), Description(description), MarketValue(0){
        }
        
        virtual ~RateInstrument(){

        }

        virtual RateInstrumentType getInstrumentType(){
            return Type;
        }

        virtual string getInstrumentDescription(){
            return Description;
        }
        
        virtual void setMarketValue(double market_value){
            MarketValue = market_value;
        }

        virtual double getMarketValue(){
            return MarketValue;
        }

        virtual double getModelValue() = 0;

        virtual double getLoss(size_t order = 1){
            double mkt_diff = getModelValue() - MarketValue;
            double loss = mkt_diff;
            size_t i;
            for ( i = 1; i < order; ++i){
                loss *= mkt_diff;
            }
            return loss;
        }

};

#endif
