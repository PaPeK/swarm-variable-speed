/*  MathTools 
    defines basic mathmatical operations on numbers and vectors 
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef mathtools_H
#define mathtools_H
#include <cmath>
#include <vector>
#include <algorithm>
#include <exception>
#include <iostream>

// definitions of macros
#define scalar(A,B) (A[0]*B[0]+A[1]*B[1])
// #define vec_length(A) (sqrt(A[0]*A[0]+A[1]*A[1]))
#define arr_length(A) (sqrt(A[0]*A[0]+A[1]*A[1]))
// #define EMPTY -1
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])

double SigmoidFunction(double x, double x0, double steepness);
void NormalizeArray(double *arr);
double *CalcDistArr(double ri[2], double rj[2], int BC, double sizeL);
// mathmatical positive rotation for vector
std::vector<double> RotateVecCcw(std::vector<double> &vec, double angle);
// mathmatical negative rotation for vector
std::vector<double> RotateVecCw(std::vector<double> &vec, double angle);
// returns intersecting postion of 2 circles
std::vector<double> CircleIntersect(double r1, double r2, double cen2x);
// returns intersecting postion of circles with line
std::vector<double> CircleLineIntersect(double r,
                                        std::vector<double> & p1,
                                        std::vector<double> & p2);
template<class T>
int sgn(T val);


template<class T>
std::vector<T> CalcDistVec(std::vector<T> &ri, std::vector<T> &rj, int BC, double sizeL)
{
    // Function for calculating the distance N-dimensional vector
    // taking into account periodic boundary conditions if BC==0
    unsigned int dim = ri.size();
    std::vector<T> r_ji(dim);
    unsigned int d;

    // Calculate connecting vector
    for(d=0; d<dim; d++)
        r_ji[d] = rj[d] - ri[d];

    if(!BC) // Check periodic boundary
        for(d=0; d<dim; d++){
            double sign = sgn(r_ji[d]);
            r_ji[d] = fmod(r_ji[d] + sign * sizeL/2, sizeL) - sign * sizeL/2;
        }
    return r_ji;
}

template<class T>
double vec_length(std::vector<T> &vec){
    unsigned int dim = vec.size();
    double len = 0;
    for (unsigned int i=0; i<dim; i++)
        len += vec[i] * vec[i];
    len = sqrt(len);
    return len;
}


template<class T>
std::vector<T> vec_sub(std::vector<T> &vec0, std::vector<T> &vec1){
    unsigned int dim = vec0.size();
    std::vector<T> vec2(dim);
    for (unsigned int i=0; i<dim; i++)
        vec2[i] = vec0[i] - vec1[i];
    return vec2;
}

// vec_add221 = adds second input vector to first input vector
template<class T, class T2>
void vec_sub221(std::vector<T> &vec1, std::vector<T2> &vec2){
    unsigned int dim = vec1.size();
    for (unsigned int i=0; i<dim; i++)
        vec1[i] -= vec2[i];
}

template<class T>
std::vector<T> vec_add(std::vector<T> &vec0, std::vector<T> &vec1){
    unsigned int dim = vec0.size();
    std::vector<T> vec2(dim);
    for (unsigned int i=0; i<dim; i++)
        vec2[i] = vec0[i] + vec1[i];
    return vec2;
}

// vec_add221 = adds second input vector to first input vector
template<class T, class T2>
void vec_add221(std::vector<T> &vec1, std::vector<T2> &vec2){
    unsigned int dim = vec1.size();
    for (unsigned int i=0; i<dim; i++)
        vec1[i] += vec2[i];
}

template<class T, class T1>
std::vector<double> vec_mul(std::vector<T> &vec, T1 factor){
    unsigned int dim = vec.size();
    std::vector<double> veco(dim);
    for (unsigned int i=0; i<dim; i++)
        veco[i] = vec[i] * factor;
    return veco;
}


template<class T>
void vec_mul221(std::vector<double> &vec, T factor){
    unsigned int dim = vec.size();
    for (unsigned int i=0; i<dim; i++)
        vec[i] *= factor;
}


template<class T>
std::vector<double> vec_set_mag(std::vector<T> &vec, double mag){
    double mag_current = vec_length(vec);
    double correction = mag / mag_current;
    std::vector<double> vec_out = vec_mul(vec, correction);
    return vec_out;
}


template<class T, class T1>
std::vector<double> vec_div(std::vector<T> &vec, T1 denom){
    unsigned int dim = vec.size();
    std::vector<double> veco(dim);
    for (unsigned int i=0; i<dim; i++)
        veco[i] = vec[i] / denom;
    return veco;
}

template<class T>
void vec_div221(std::vector<double> &vec, T denom){
    unsigned int dim = vec.size();
    for (unsigned int i=0; i<dim; i++)
        vec[i] /= denom;
}

template<class T>
T vec_sum(std::vector<T> &vec){
    unsigned int dim = vec.size();
    T sum = 0;
    for (unsigned int i=0; i<dim; i++)
        sum += vec[i];
    return sum;
}

template<class T, class T2>
double vec_dot(std::vector<T> &vec, std::vector<T2> &vec2){
    unsigned int dim = vec.size();
    double dot = 0;
    for (unsigned int i=0; i<dim; i++)
        dot += vec[i] * vec2[i];
    return dot;
}


template<class T>
std::vector<T> vec_perp(std::vector<T> &vec){
    std::vector<T> veco(2);
    veco[0] = vec[1];
    veco[1] = - vec[0];
    return veco;
}


template<class T>
double CalcDist(std::vector<T> &ri, std::vector<T> &rj, int BC, double sizeL)
{
    // Function for calculating the distance N-dimensional vector
    // taking into account periodic boundary conditions if BC==0
    unsigned int dim = ri.size();
    std::vector<T> r_ji(dim);
    unsigned int d;

    // Calculate connecting vector
    for(d=0; d<dim; d++)
        r_ji[d] = rj[d] - ri[d];

    if(!BC) // Check periodic boundary
        for(d=0; d<dim; d++){
            double sign = sgn(r_ji[d]);
            r_ji[d] = fmod(r_ji[d] + sign * sizeL/2, sizeL) - sign * sizeL/2;
        }
    double dist = 0;
    for(d=0; d<dim; d++)
        dist += r_ji[d] * r_ji[d];
    return sqrt(dist);
}


template<class T>
int where_val_in_vector(std::vector<T> &vec, T a){
        // Function returns index where the integer a is first located in vector vec
        // if a is not found, returns size of array
        // std::vector<T>::iterator it;
        // it = std::find(vec.begin(), vec.end(), a);
        // int dist = std::distance(vec.begin(), it);
        int dist = std::distance(vec.begin(), std::find(vec.begin(), vec.end(), a));
        return dist;
}


// returns -1, 0, 1 accroding to signum
template<class T>
int sgn(T val){
    return (T(0) < val) - (val < T(0));
}


template<class T>
double Get_quantile(std::vector<T> &data, double quantile){
    unsigned int N = data.size();
    unsigned int q = quantile * N;
    std::sort(data.begin(), data.end());
    return data[q];
}

// Comparison logic needed for partial-sort applied in FOV-code
template <typename T>
struct Comp{
            Comp( const std::vector<T>& v ) : _v(v) {}
                bool operator ()(T a, T b) { return _v[a] > _v[b]; }
                    const std::vector<T>& _v;
};

template <class T, class T1>
double DistPoint2Segment(std::vector<T> &x, std::vector<T1> &A, std::vector<T1> &B){
    std::vector<double> Ax {x[0] - A[0], x[1] - A[1]};
    std::vector<double> Bx {x[0] - B[0], x[1] - B[1]};
    std::vector<double> AB {B[0] - A[0], B[1] - A[1]};
    double dAx = vec_length(Ax);
    double dBx = vec_length(Bx);
    double dist;
    double AB_Ax = vec_dot(AB, Ax);
    double AB_Bx = vec_dot(AB, Bx);
    if (dAx < 0.0001 || dBx < 0.0001){
        dist = 0;
    }
    // case: dot(AB, Ax)<0
    //    x
    //     \
    // -----A-----
    //      |
    //      |
    // -----B-----
    else if (AB_Ax < 0)
        dist = dAx;
    // case: dot(AB, Bx)>0
    // -----A-----
    //      |
    //      |
    // -----B-----
    //     /
    //    x
    else if (AB_Bx > 0)
        dist = dBx;
    // case: dot(AB, Ax)>0 and dot(AB, Bx)<0
    //          AND AB  is vertical 
    // -----A-----
    //      |
    //  x---|
    // -----B-----
    else if (AB[0] == 0){
        dist = abs(Ax[0]);
        // if (abs(Ax[0] - Bx[0]) > 0.0001)
        //     throw std::exception::runtime_error;
    }
    // case: dot(AB, Ax)>0 and dot(AB, Bx)<0
    //          AND AB is horizontal 
    // -----A-----B-----
    //          |
    //          |
    //          x
    else if (AB[1] == 0){
        dist = abs(Ax[1]);
        // if (abs(Ax[1] - Bx[1]) > 0.0001)
        //     throw std::runtime_error;
    }
    // case: dot(AB, Ax)>0 and dot(AB, Bx)<0
    //          AND not AB is not vertical or horizontal
    // --------A-
    //        /
    //  x----/
    // -----B----
    else{
        // formulating line as : ax + by + c = 0
        double a = A[1] - B[1];
        double b = B[0] - A[0];
        double c = A[0] * B[1] - B[0] * A[1];
        dist = abs(a * x[0] + b * x[1] + c) / sqrt(a * a + b * b);
    }
    return dist;
}

#endif
