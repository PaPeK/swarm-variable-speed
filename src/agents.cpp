#include "agents.h"


std::vector<double> particle::out(void){
    std::vector<double> out;
    out.reserve(7);
    out.push_back( x[0]); 
    out.push_back( x[1]); 
    out.push_back( v[0]); 
    out.push_back( v[1]); 
    out.push_back( fitness); 
    out.push_back( force[0]); 
    out.push_back( force[1]); 
    // out.push_back( force_flee[0]); 
    // out.push_back( force_flee[1]); 
    return out;
}


std::vector<double> predator::out(void){
    std::vector<double> out;
    out.reserve(4);
    out.push_back( x[0]); 
    out.push_back( x[1]); 
    out.push_back( v[0]); 
    out.push_back( v[1]); 
    return out;
}
