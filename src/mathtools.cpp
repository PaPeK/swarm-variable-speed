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
#include "mathtools.h"
double SigmoidFunction(double x, double x0, double steepness)
{
    // Simple sigmoid function going from x->[0,1]
    // x0 sets the transition point f(x0)=0.5;
    // steepness sets the sharpness of the transition
    // steepness > 0: transition from 0 to 1
    // steepness < 0: transition from 1 to 0

    double y;
    y=0.5*(tanh(steepness*(x-x0))+1.0);
    return y;
}

void NormalizeArray(double *arr)
{
    // Normalize the input vector to unit length
    double vnorm = arr_length(arr);
    if(vnorm>0)
    {
        arr[0]/=vnorm;
        arr[1]/=vnorm;
    }
}

double *CalcDistArr(double ri[2], double rj[2], int BC, double sizeL)
{
    // Function for calculating the distance vector
    // taking into account periodic boundary conditions if BC==0
    static double r_ji[2];
    int d;

    // Calculate connecting vector
    for(d=0;d<2;d++)
        r_ji[d]= 
                rj[d]-
                ri[d];

    if(!BC) // Check periodic boundary
        for(d=0;d<2;d++)
            r_ji[d]=fmod(r_ji[d] + sizeL/2,sizeL)-sizeL/2;
    return r_ji;
}

std::vector<double> RotateVecCcw(std::vector<double> &vec, double angle){
    // rotation of vec in 2D in mathmatical positive direction (ccw)
    std::vector<double> vec_r(2);
    vec_r[0] = cos(angle) * vec[0] - sin(angle) * vec[1];
    vec_r[1] = sin(angle) * vec[0] + cos(angle) * vec[1];
    return vec_r;
}

std::vector<double> RotateVecCw(std::vector<double> &vec, double angle){
    // rotation of vec in 2D in mathmatical negative direction (cw)
    std::vector<double> vec_r(2);
    vec_r[0] = cos(angle) * vec[0] + sin(angle) * vec[1];
    vec_r[1] = - sin(angle) * vec[0] + cos(angle) * vec[1];
    return vec_r;
}


std::vector<double> CircleIntersect(double r1, double r2, double cen2x){
    // assumes that:
    // - center of circle 1 is (0, 0)
    // - center of circle 2 is on y axis: (cen2x, 0)
    std::vector<double> inter(2);
    // from using both circle equations: x^2 + y^2 = R^2, and (x-cen2x)^2 + y^2 = r^2
    inter[0] = ( r1*r1 - r2*r2 + cen2x*cen2x ) / ( 2 * cen2x );
    //pluggin in into x^2 + y^2 = R^2 results in:
    inter[1] = sqrt(r1*r1 - inter[0]*inter[0]);
    return inter;
}


std::vector<double> CircleLineIntersect(double r,
                                        std::vector<double> & p1,
                                        std::vector<double> & p2){
    // equations from http://mathworld.wolfram.com/Circle-LineIntersection.html
    // assumes that:
    // - center of circle at origin (0, 0)
    // - p1 and p2 are enpoints of line-segment and desired 
    //      intersection lies in between there 2
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dr = sqrt( dx * dx + dy * dy );
    double D = p1[0] * p2[1] - p2[0] * p1[1];
    // define borders to only return 1 intersection which lies in between p1 and p1
    double lx = fmin(p1[0], p2[0]);
    double ux = fmax(p1[0], p2[0]);
    double ly = fmin(p1[1], p2[1]);
    double uy = fmax(p1[1], p2[1]);
    double x0, y0;
    x0 = ( D * dy + sgn(dy) * dx * sqrt( r*r * dr*dr - D*D ) ) / (dr*dr);
    if ( (lx < x0) && (x0 < ux) )
        y0 = ( -D * dx + fabs(dy) * sqrt( r*r * dr*dr - D*D ) ) / (dr*dr);
    else{
        x0 = ( D * dy - sgn(dy) * dx * sqrt( r*r * dr*dr - D*D ) ) / (dr*dr);
        y0 = ( -D * dx - fabs(dy) * sqrt( r*r * dr*dr - D*D ) ) / (dr*dr);
    }
    // // debugging:
    // if (x0 < lx || ux < x0 || y0 < ly || uy < y0){
    //     std::cout<< "ERROR intersection not between p1 and p2" << std::endl;
    //     std::cout<< "mi " << lx << " " << ly <<std::endl;
    //     std::cout<< "xy " << x0 << " " << y0 <<std::endl;
    //     std::cout<< "ma " << ux << " " << uy <<std::endl;
    // }
    std::vector<double> out {x0, y0};
    return out;
}
