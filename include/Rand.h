/*
 * ==========================================================================
 *
 *       Filename:  Rand.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-09-16
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jimmy1500
 *                  lighteningmagic@gmail.com
 *
 * ==========================================================================
 */

#ifndef RAND_H
#define RAND_H

//#include <iostream>
//#include <iomanip>
//#include <string>
//#include <map>
#include <random>

using namespace std;
//This class creates a 2D pool of random number engines with ndims & npaths <user definited> for multi-threaded use
//Warning: do not share a particular random engine between objects in multithreaded settings

#ifdef WIN32 
    typedef mt19937 Engine;
#else
    typedef mt19937_64 Engine;
#endif 
        
template<class T>
class Rand{
    private:
        size_t ndims;
        size_t npaths;
        Engine ** Engines;
    public:

        normal_distribution<T> norm_dist;
        uniform_real_distribution<T> uniform_dist;
        
        //Constructor
        Rand(size_t, size_t, T, T); //instantiate object with random_device produced seeds (system entropy)
        Rand(size_t, size_t, T, T, unsigned int**);//instantiate object with user provided seeds
        
        //Destructor
        ~Rand(void);

        //Resize generator pool and re-seed new random engines using system entropy
        size_t resize(size_t, size_t);

        //generate uniform random i.i.d using (ndim, npath)th generator in the pool 
        T urand(size_t, size_t);
        //generate normal random i.i.d using (ndim, npath)th generator in the pool 
        T nrand(size_t, size_t);
        
        //utilities
        void seed();//seed/re-seed all generators in the pool with random_device produced intergers
        void seed(unsigned int**);//seed/re-seed all generators with user provided interger seeds
        
        Engine& getEngine(size_t, size_t);
        Engine* operator[](size_t);

        size_t getPaths();
        size_t getDimensions();
};


template <class T>
Rand<T>::Rand(size_t paths, size_t dims, T mu, T sigma) : ndims(dims), npaths(paths), Engines(new Engine* [dims]), norm_dist(mu, sigma), uniform_dist(mu, sigma)
{
    random_device rd{};

    size_t i,j;
    for (i = 0; i < ndims; ++i){
        Engines[i] = new Engine[npaths];
        for (j = 0; j < npaths; ++j){
            Engines[i][j] = Engine(rd()); //seed engine with system entropy
        }
    } 
}

template<class T>
Rand<T>::Rand(size_t paths, size_t dims, T mu, T sigma, unsigned int**rseed) : Engines(new Engine* [dims]), npaths(paths), ndims(dims), norm_dist(mu, sigma), uniform_dist(mu, sigma)
{
    size_t i,j;
    for (i = 0; i < ndims; ++i){
        Engines[i] = new Engine[npaths];
        for (j = 0; j < npaths; ++j){
            Engines[i][j] = Engine(rseed[i][j]); //seed engine with user defined seeds
        }
    }
}

template <class T>
Rand<T>::~Rand(void){
    size_t i;
    for (i = 0; i < ndims; ++i){
        delete [] Engines[i];
    }
    delete [] Engines;
}

template <class T>
size_t Rand<T>::resize(size_t paths, size_t dims){
    if ( ndims != dims ){
        size_t i,j;
        for (i = 0; i < ndims; ++i){
            delete [] Engines[i];
        }
        delete [] Engines; Engines = new Engine* [dims];
        random_device rd{};
        for (i = 0; i < dims; ++i){
            Engines[i] = new Engine[paths];
            for (j = 0; j < paths; ++j){
                Engines[i][j] = Engine(rd()); //seed engine with system entropy
            }
        }
        npaths = paths; ndims = dims;
        return 1;
    }else if ( npaths != paths ){
        size_t i,j;
        random_device rd{};
        for (i = 0; i < ndims; ++i){
            delete [] Engines[i]; Engines[i] = new Engine[paths];
            for (j = 0; j < paths; ++j){
                Engines[i][j] = Engine(rd()); //seed engine with system entropy
            }
        }
        npaths = paths;
        return 1;
    }
    return 0;
}

template<class T>
T Rand<T>::urand(size_t dim, size_t path){
    return uniform_dist(Engines[dim][path]);
}

template<class T>
T Rand<T>::nrand(size_t dim, size_t path){//Access Random Generators Pool by id (id<=engine_pool_size), return norm dist value
   return norm_dist(Engines[dim][path]);
}

template<class T>
void Rand<T>::seed(){
    random_device rd{};

    size_t i,j;
    for (i = 0; i < ndims; ++i){
        for (j = 0; j < npaths; ++j){
            Engines[i][j].seed(rd());
        }
    }
}

template<class T>
void Rand<T>::seed(unsigned int**rseed){
    size_t i,j;
    for (i = 0; i < ndims; ++i){
        for (j = 0; j < npaths; ++j){
            Engines[i][j].seed(rseed[i][j]);
        }
    }
}

template<class T>
Engine& Rand<T>::getEngine(size_t dim, size_t path){
    return Engines[dim][path];
}

template<class T>
Engine* Rand<T>::operator[](size_t dim){
    return Engines[dim];
}

template<class T>
size_t Rand<T>::getPaths(){
    return npaths;
}

template<class T>
size_t Rand<T>::getDimensions(){
    return ndims;
}

#endif //RAND_H
