/*
 * ==========================================================================
 *
 *       Filename:  Matrix.h
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
#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <algorithm>

using namespace std;

struct mat3d{
    size_t nterms;
    size_t npaths;
    size_t ndims;
    double *** value;

    mat3d(mat3d & mtrx){
        nterms = mtrx.nterms;
        npaths = mtrx.npaths;
        ndims  = mtrx.ndims;
        value  = mtrx.value;
    }

    //allocate and zero memory dynamically
    mat3d(size_t terms, size_t paths, size_t dims)
        : nterms(terms),
        npaths(paths),
        ndims(dims),
        value(new double **[dims]){

        size_t dim, path;
        for (dim=0; dim<ndims; ++dim){
            value[dim] = new double*[npaths];
            for (path=0; path<npaths; ++path){
                value[dim][path] = new double [nterms](); //zero memory
            }
        }
    }
    
    //delete matrix memory
    ~mat3d(){
        size_t dim, path;
        for (dim=0; dim<ndims; ++dim){
            for (path=0; path<npaths; ++path){
                delete [] value[dim][path];
            }
            delete [] value[dim];
        }
        delete [] value;
    }
    
    //operator overload
    double ** operator[](size_t dim){
        return value[dim];
    }
    
    size_t resize(size_t terms, size_t paths, size_t dims){
        if (ndims != dims){
            size_t dim, path;
            for (dim=0; dim<ndims; ++dim){
                for (path=0; path<npaths; ++path){
                    delete [] value[dim][path];
                }
                delete [] value[dim];
            }
            delete [] value; value = new double **[dims];
            for (dim=0; dim<dims; ++dim){
                value[dim] = new double*[paths];
                for (path=0; path<paths; ++path){
                    value[dim][path] = new double [terms](); //zero memory
                }
            }
            nterms = terms; npaths = paths; ndims = dims;
            return 1;
        }else if (npaths != paths){
            size_t dim, path;
            for (dim=0; dim<ndims; ++dim){
                for (path=0; path<npaths; ++path){
                    delete [] value[dim][path];
                }
                delete [] value[dim]; value[dim] = new double*[paths];
                for (path=0; path<paths; ++path){
                    value[dim][path] = new double [terms](); //zero memory
                }
            }
            nterms = terms; npaths = paths;
            return 1;
        }else if (nterms != terms){
            size_t dim, path;
            for (dim=0; dim<ndims; ++dim){
                for (path=0; path<npaths; ++path){
                    delete [] value[dim][path]; value[dim][path] = new double [terms](); //zero memory
                }
            }
            nterms = terms;
            return 1;
        }
        return 0;
    }
};

struct mat2d{
    size_t nterms;
    size_t npaths;
    double ** value;
    
    mat2d(mat2d & mtrx){
        nterms = mtrx.nterms;
        npaths = mtrx.npaths;
        value  = mtrx.value;
    }

    //allocate and zero memory dynamically
    mat2d(size_t terms, size_t paths)
        : nterms(terms),
        npaths(paths),
        value(new double *[paths]){

        size_t path;
        for (path=0; path<npaths; ++path){
            value[path] = new double [nterms](); //zero memory
        }
    }
    
    //delete matrix memory
    ~mat2d(){
        size_t path;
        for (path=0; path<npaths; ++path){
            delete [] value[path];
        }
        delete [] value;
    }
    
    //operator overload
    double * operator[](size_t path){
        return value[path];
    }

    size_t resize(size_t terms, size_t paths){
        if (npaths != paths){
            size_t path;
            for (path=0; path<npaths; ++path){
                delete [] value[path];
            }
            delete [] value; value = new double *[paths];
            for (path=0; path<paths; ++path){
                value[path] = new double [terms](); //zero memory
            }
            nterms = terms; npaths = paths;
            return 1;
        } else if (nterms != terms){
            size_t path;
            for (path=0; path<npaths; ++path){
                delete [] value[path]; value[path] = new double [terms](); //zero memory
            }
            nterms = terms;
            return 1;
        }
        return 0;
    }
}; 

struct mat1d{
    size_t length;
    double * value;

    mat1d(mat1d & vector){
        length = vector.length;
        value  = vector.value;
    }

    mat1d(size_t len) : length(len), value(new double [len]()){}

    ~mat1d(){
        delete [] value;
    }

    //operator overload
    double operator[](size_t idx){
        return value[idx];
    }

    size_t resize(size_t len, bool refresh = true){
        if (len != length){
            delete [] value; value = new double[len]();
            length = len;
            return 1;
        }else if (refresh){ //only necessary when value is not instantiated
            fill(value, value+length, 0);
            return 1;
        }
        return 0;
    }
};

class Matrix{
    mat2d mat;
    Matrix(size_t n_row, size_t n_col) : mat(n_row, n_col) {

    }

    size_t nColns(){
        return mat.npaths;
    }

    size_t nRows(){
        return  mat.nterms;
    }

    double * operator [] (size_t col){
        return mat[col];
    }

    double operator() (size_t row, size_t col){
        return mat[col][row];
    }

    Matrix operator* (Matrix & B){
        if (this->nColns() == B.nRows()){
            Matrix C = Matrix(this->nRows(), B.nColns());

            size_t col, row, i;
            for ( col = 0; col < C.nColns(); ++col ){
                for ( row = 0; row < C.nRows(); ++row ){
                    for ( i = 0; i < B.nRows(); ++i){
                        C[col][row] += B[col][i]*(*this[i][row]);
                    }
                }
            }
            return C;
        }else{
            return *this;
        }
    }

    Matrix Inverse(){
        return *this;
    }
};

#endif
