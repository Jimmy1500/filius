#include "Curve.h"

Curve::Curve(Curve * curve){ //Shallow
    *this = &(*curve);
}

Curve::Curve(double * terms, double * prices, size_t len) //Deep
    : Terms(new double[len]),
    Values(new double[len]),
    Length(len){
    DeepCopy_YC(terms, prices)
}

Curve::~Curve(){
    delete [] Terms;
    delete [] Values;
}

void Curve::interpolateTerms(double t, size_t& index, int& direction){
    if ( !Length ){ throw 0; }

    size_t lower = 0, upper = Length-1;
    if ( upper == lower || t == Terms[lower] ){ index = lower; direction = 0; return; }
    if ( t < Terms[lower] ){ index = lower; direction = 1; return; }
    if ( t == Terms[upper] ){ index = upper; direction = 0; return; }
    if ( t > Terms[upper] ){ index = upper; direction = -1; return; }

    while ( upper-lower > 1 ){
        index = ( lower+upper ) / 2;
        if ( t == Terms[index] ){ direction = 0; return; }
        else if ( t < Terms[index] ){ upper = index; }
        else { lower = index; }
    }

    index = lower; direction = 1; return;
}

double Curve::P(double tau){
    if (tau < 0.0){ tau = -tau; } //Assume past mirrors current, do not abuse!!
    else if (tau == 0.0){ return 1.0; }

    size_t t, T; int direction; interpolateTerms(tau, t, direction);
    if (!direction){ return Values[t]; } else{ T = t + direction; }

    // log linear: lnP(time, Tau) = ( (T-Tau)*lnP(time,t) + (Tau-t)*lnP(time,T) ) / (T-t)
    return exp( ((Terms[T]-tau) * log(Values[t]) +
                (tau-Terms[t]) * log(Values[T])) / (Terms[T]-Terms[t]) );
}

double Curve::P(double t, double T){
    if (t == T){
        return 1.0;
    }else{
        return P(T)/P(t);
    }
}

double * Curve::getTerms(){
    return Terms;
}

double * Curve::getValues(){
    return Values;
}

size_t Curve::getLength(){
    return Length;
}
