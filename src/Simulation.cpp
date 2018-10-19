#include "Simulation.h"

Simulation::Simulation():mu(0.0),sigma(1.0),Generator(nullptr){
}

Simulation::Simulation(double mu, double sigma):mu(mu),sigma(sigma),Generator(nullptr){
}

Simulation::~Simulation(){
    delete Generator;
}

//private
//make sure generator has the same dimensions as those of the random matrix specified by user
void Simulation::fitGenerator(size_t npaths, size_t ndims){
    if (Generator){//check if generator has been instantiated
        Generator->resize(npaths, ndims);
    }else{
        Generator = new Rand<double>(npaths, ndims, mu, sigma);
    }
}

//public
//--------
void Simulation::setRandAttributes(double mu, double sigma){
   this->mu = mu;
   this->mu = sigma;
}

//Multithreaded random matrix generation processes:
//thread needs an instance of Simulation to operate member function on, here we used lambda expression
//lambda expression:
//[a,&b] where a is captured by value and b is captured by reference.
//[this] captures the this pointer by value
//[&] captures all automatic variables odr-used in the body of the lambda by reference
//[=] captures all automatic variables odr-used in the body of the lambda by value
//[] captures nothing/
void Simulation::ugen(mat3d & rand_mtrx, size_t nthreads, bool match_dims){
    if (nthreads < 1){
       nthreads = thread::hardware_concurrency();
    }
    if (nthreads > rand_mtrx.npaths){
       nthreads = rand_mtrx.npaths;
    }

    switch (nthreads){
        case 1: 
        // optimized for single threade(use main thread to avoid overhead)
           if (match_dims){
              fitGenerator(rand_mtrx.npaths, rand_mtrx.ndims);
              size_t dim, path, term;
              for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                 for (path = 0; path < rand_mtrx.npaths; ++path){
                    for (term = 0; term < rand_mtrx.nterms; ++term){
                       rand_mtrx[dim][path][term] = Generator->urand(dim, path);
                    }
                 }
              }
           }else{
              fitGenerator(1, 1);
              size_t dim, path, term;
              for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                 for (path = 0; path < rand_mtrx.npaths; ++path){
                    for (term = 0; term < rand_mtrx.nterms; ++term){
                       rand_mtrx[dim][path][term] = Generator->urand(0, 0);
                    }
                 }
              }
           }
           break;

        default:
            size_t extra_paths = rand_mtrx.npaths % nthreads;
            size_t paths_per_thread = (rand_mtrx.npaths - extra_paths)/nthreads;
            size_t last_thread = nthreads-1;
            size_t t;

            thread * threads = new thread[nthreads];
            if (match_dims){
               fitGenerator(rand_mtrx.npaths, rand_mtrx.ndims);
               for (t = 0; t < last_thread; ++t){
                  threads[t] = thread([this,t,paths_per_thread,&rand_mtrx]()->void{
                        size_t begin = t*paths_per_thread;
                        size_t end = (t+1)*paths_per_thread;
                        size_t dim, path, term;
                        for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                           for (path = begin; path < end; ++path){
                              for (term = 0; term < rand_mtrx.nterms; ++term){
                                 rand_mtrx[dim][path][term] = Generator->urand(dim, path);
                              }
                           }
                        }
                  });
               }

               threads[last_thread] = thread([this,last_thread,paths_per_thread,&rand_mtrx]()->void{
                     size_t begin = last_thread*paths_per_thread;
                     size_t dim, path, term;
                     for (dim=0; dim < rand_mtrx.ndims; ++dim){
                        for (path = begin; path < rand_mtrx.npaths; ++path){
                           for (term = 0; term < rand_mtrx.nterms; ++term){
                              rand_mtrx[dim][path][term] = Generator->urand(dim, path);
                           }
                        }
                     }
               });

            }else{
               fitGenerator(nthreads, 1);
               for (t = 0; t < last_thread; ++t){
                  threads[t] = thread([this,t,paths_per_thread,&rand_mtrx]()->void{
                        size_t begin = t*paths_per_thread;
                        size_t end = (t+1)*paths_per_thread;
                        size_t dim, path, term;
                        for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                           for (path = begin; path < end; ++path){
                              for (term = 0; term < rand_mtrx.nterms; ++term){
                                 rand_mtrx[dim][path][term] = Generator->urand(0, t);
                              }
                           }
                        }
                  });
               }

               threads[last_thread] = thread([this,last_thread,paths_per_thread,&rand_mtrx]()->void{
                     size_t begin = last_thread*paths_per_thread;
                     size_t dim, path, term;
                     for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                        for (path = begin; path < rand_mtrx.npaths; ++path){
                           for (term = 0; term < rand_mtrx.nterms; ++term){
                              rand_mtrx[dim][path][term] = Generator->urand(0, last_thread);
                           }
                        }
                     }
               });
            }
            
            for (t=0; t < nthreads; ++t){
                threads[t].join();
            }
            
            thread garbage_collection = thread([threads]()->void{ delete [] threads; }); 
            garbage_collection.detach();
            break;
    }
}

void Simulation::ngen(mat3d & rand_mtrx, size_t nthreads, bool match_dims){//nthreads<=rand_mtrx.npaths
    if (nthreads < 1){
       nthreads = thread::hardware_concurrency();
    }
    if (nthreads > rand_mtrx.npaths){
       nthreads = rand_mtrx.npaths;
    }

    switch (nthreads){
        case 1: 
        // optimized for single threade(use main thread to avoid overhead)
           if (match_dims){
              fitGenerator(rand_mtrx.npaths, rand_mtrx.ndims);
              size_t dim, path, term;
              for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                 for (path = 0; path < rand_mtrx.npaths; ++path){
                    for (term = 0; term < rand_mtrx.nterms; ++term){
                       rand_mtrx[dim][path][term] = Generator->nrand(dim, path);
                    }
                 }
              }
           }else{
              fitGenerator(1, 1);
              size_t dim, path, term;
              for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                 for (path = 0; path < rand_mtrx.npaths; ++path){
                    for (term = 0; term < rand_mtrx.nterms; ++term){
                       rand_mtrx[dim][path][term] = Generator->nrand(0, 0);
                    }
                 }
              }
           }
           break;

        default:
            size_t extra_paths = rand_mtrx.npaths % nthreads;
            size_t paths_per_thread = (rand_mtrx.npaths - extra_paths)/nthreads;
            size_t last_thread = nthreads-1;
            size_t t;

            thread * threads = new thread[nthreads];
            if (match_dims){
               fitGenerator(rand_mtrx.npaths, rand_mtrx.ndims);
               for (t = 0; t < last_thread; ++t){
                  threads[t] = thread([this,t,paths_per_thread,&rand_mtrx]()->void{
                        size_t begin = t*paths_per_thread;
                        size_t end = (t+1)*paths_per_thread;
                        size_t dim, path, term;
                        for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                           for (path = begin; path < end; ++path){
                              for (term = 0; term < rand_mtrx.nterms; ++term){
                                 rand_mtrx[dim][path][term] = Generator->nrand(dim, path);
                              }
                           }
                        }
                  });
               }

               threads[last_thread] = thread([this,last_thread,paths_per_thread,&rand_mtrx]()->void{
                     size_t begin = last_thread*paths_per_thread;
                     size_t dim, path, term;
                     for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                        for (path = begin; path < rand_mtrx.npaths; ++path){
                           for (term = 0; term < rand_mtrx.nterms; ++term){
                              rand_mtrx[dim][path][term] = Generator->nrand(dim, path);
                           }
                        }
                     }
               });

            }else{
               fitGenerator(nthreads, 1);
               for (t = 0; t < last_thread; ++t){
                  threads[t] = thread([this,t,paths_per_thread,&rand_mtrx]()->void{
                        size_t begin = t*paths_per_thread;
                        size_t end = (t+1)*paths_per_thread;
                        size_t dim, path, term;
                        for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                           for (path = begin; path < end; ++path){
                              for (term = 0; term < rand_mtrx.nterms; ++term){
                                 rand_mtrx[dim][path][term] = Generator->nrand(0, t);
                              }
                           }
                        }
                  });
               }

               threads[last_thread] = thread([this,last_thread,paths_per_thread,&rand_mtrx]()->void{
                     size_t begin = last_thread*paths_per_thread;
                     size_t dim, path, term;
                     for (dim = 0; dim < rand_mtrx.ndims; ++dim){
                        for (path = begin; path < rand_mtrx.npaths; ++path){
                           for (term = 0; term < rand_mtrx.nterms; ++term){
                              rand_mtrx[dim][path][term] = Generator->nrand(0, last_thread);
                           }
                        }
                     }
               });
            }
            
            for (t = 0; t < nthreads; ++t){
                threads[t].join();
            }
            
            thread garbage_collection = thread([threads]()->void{ delete [] threads; }); 
            garbage_collection.detach();
            break;
    }
}

Rand<double> * Simulation::getGenerator(){
    return Generator;
}
