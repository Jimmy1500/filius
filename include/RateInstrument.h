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
 *         Author:  Jimmy1500
 *                  lighteningmagic@gmail.com
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

        virtual double getMarketSpread(size_t order = 1){
            double spread = getModelValue() - MarketValue;
            size_t i;
            for ( i = 1; i < order; ++i){
                spread *= spread;
            }
            return spread;
        }

};

#endif
