#include "agents.h"

std::vector<double> particle::out(void){
    std::vector<double> out;
    out.reserve(6);
    out.push_back( x[0]); 
    out.push_back( x[1]); 
    out.push_back( v[0]); 
    out.push_back( v[1]); 
    out.push_back( force[0]); 
    out.push_back( force[1]); 
    return out;
}
