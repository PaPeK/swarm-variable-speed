/*  AgentsDynamics
    defines movement rules/differential equations of agents in
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser
*/ 

#include "agents_dynamics.h"


void MoveParticle(particle &a, params * ptrSP, gsl_rng *r, double rnp){
    double h1 = 0;      // helper
    double dt = ptrSP->dt;
    double rep_strength = ptrSP->rep_strength;
    double forcep = 0.0;
    double forcev = 0.0;
    std::vector<double> force(2, 0);
    std::vector<double> hvec(2, 0);
    double lphi = 0.0;
    int BC = ptrSP->BC;
    double sizeL = ptrSP->sizeL;
    double alg_strength = ptrSP->alg_strength;
    double turn_alpha = ptrSP->turn_alpha;
    double beta = ptrSP->beta;
    double speed0 = ptrSP->speed0;

    // Calc total social force
    if( a.counter_rep > 0 ){
        vec_mul221(a.force_rep, rep_strength / a.counter_rep);
        vec_mul221(a.force_alg, alg_strength / a.counter_rep);
        vec_add221(force, a.force_rep);
        vec_add221(force, a.force_alg);
    }
    a.force = force; // only for output
    
    // Calculate polar angle
    lphi = a.phi;
    double vproj = a.vproj;     // to use correct time-step
    forcev = force[0] * cos(lphi) + force[1] * sin(lphi);
    a.vproj += (beta * (speed0 - a.vproj) + forcev) * dt;
    // a.vproj += rnv;     // rnv = sqrt(dt * Dv) * N(0, 1) TODO: should we have speed noise?
    if (a.vproj < 0)     // prevents F of swimming back
        a.vproj = 0;    // angle adapted below to exactly the force direction

    forcep =- force[0] * sin(lphi) + force[1] * cos(lphi);

    if (a.vproj != 0) // TODO: should be True also for small vproj (otherwise angle overshoot) -> what is small (depends on dt and force strength)
        lphi += ( forcep * dt + rnp) / (a.vproj + turn_alpha);  // rnp = sqrt(dt * Dphi) * N(0, 1) (Wiener Process) 
    else
        lphi = atan2(force[1], force[0]); // instantaneous direction adaption
    lphi = fmod(lphi, 2*M_PI);
    lphi = fmod(2*M_PI + lphi, 2*M_PI);   // to ensure positive angle definition
    a.phi = lphi;
    a.u[0] = cos(lphi);
    a.u[1] = sin(lphi);
    
    // Move particles with speed in units of [vel.al. range. / time]
    a.v = vec_mul(a.u, a.vproj);
    a.x[0] += a.v[0]*dt;
    a.x[1] += a.v[1]*dt;
    
    // Reset all forces
    // a.force[0]    =a.force[1]=0.0;
    // for output: use effective force:
    a.force_rep[0] = a.force_rep[1] = 0.0;
    a.force_alg[0] = a.force_alg[1] = 0.0;
    a.counter_rep = 0;
    
    // Take into account boundary
    if(BC!=-1)
        Boundary(a, sizeL, BC);
}


template<class agent>
void Boundary(agent &a, double sizeL,  int BC)
{
    // Function for calculating boundary conditions
    // and update agent position and velocity accordingly
    double tx;
    double dphi=0.1;
    double tmpphi;
    double dx=0.0;
    double dist2cen;
    double diff;
    std::vector<double> hv(2);
    std::vector<double> wall_normal = a.x;
    // -1 is Open boundary condition
    switch (BC) {
        case 0:
            // Periodic boundary condition
            a.x[0] = fmod(a.x[0] + sizeL, sizeL);
            a.x[1] = fmod(a.x[1] + sizeL, sizeL);
            break;
        // TODO: any case which changes the velocity should also change u, phi
        case 1:
            // Inelastic box boundary condition
            if (a.x[0]>sizeL)
            {
                a.x[0]=0.9999*sizeL;
                if(a.v[1]>0.)
                    tmpphi=(0.5+dphi)*M_PI;
                else
                    tmpphi=-(0.5+dphi)*M_PI;
    
                a.v[0]=cos(tmpphi);
                a.v[1]=sin(tmpphi);
            }
            else if (a.x[0]<0)
            {
                a.x[0]=0.0001;
                if(a.v[1]>0.)
                    tmpphi=(0.5-dphi)*M_PI;
                else
                    tmpphi=-(0.5-dphi)*M_PI;
    
    
                a.v[0]=cos(tmpphi);
                a.v[1]=sin(tmpphi);
            }
            if (a.x[1]>sizeL)
            {
                a.x[1]=0.9999*sizeL;
                if(a.v[0]>0.)
                    tmpphi=-dphi;
                else
                    tmpphi=M_PI+dphi;
    
                a.v[0]=cos(tmpphi);
                a.v[1]=sin(tmpphi);
    
            }
            else if (a.x[1]<0)
            {
                a.x[1]=0.0001;
                if(a.v[0]>0.)
                    tmpphi=dphi;
                else
                    tmpphi=M_PI-dphi;
    
                a.v[0]=cos(tmpphi);
                a.v[1]=sin(tmpphi);
            }
    
    
            break;
        case 2:
            // Elastic box boundary condition
            if (a.x[0]>sizeL)
            {
                dx=2.*(a.x[0]-sizeL);
                a.x[0]-=dx;
                a.v[0]*=-1.;
            }
            else if (a.x[0]<0)
            {
                dx=2.*a.x[0];
                a.x[0]-=dx;
                a.v[0]*=-1.;
            }
            if (a.x[1]>sizeL)
            {
                dx=2.*(a.x[1]-sizeL);
                a.x[1]-=dx;
                a.v[1]*=-1.;
            }
            else if (a.x[1]<0)
            {
                dx=2.*a.x[1];
                a.x[1]-=dx;
                a.v[1]*=-1.;
            }
            break;
    
        case 3:
            // Periodic boundary condition in x
            // Elastic  boundary condition in y
            tx = fmod(a.x[0], sizeL);
            if (tx < 0.0f)
                tx += sizeL;
            a.x[0]=tx;
    
            if (a.x[1]>sizeL)
            {
                dx=2.*(a.x[1]-sizeL);
                a.x[1]-=dx;
                a.v[1]*=-1.;
            }
            else if (a.x[1]<0)
            {
                dx=2.*a.x[1];
                a.x[1]-=dx;
                a.v[1]*=-1.;
            }
            break;
        case 4:
            // Periodic boundary condition in x
            // Inelastic  boundary condition in y
            tx = fmod(a.x[0], sizeL);
            if (tx < 0.0f)
            {
                tx += sizeL;
            }
            a.x[0]=tx;
    
            if (a.x[1]>sizeL)
            {
                dx=(a.x[1]-sizeL);
                a.x[1]-=dx+0.0001;
                a.v[1]=0.0;
            }
            else if (a.x[1]<0)
            {
                dx=a.x[1];
                a.x[1]-=dx-0.0001;
                a.v[1]=0.0;
            }
            break;
        case 5:
            // Elastic circle boundary condition
            dist2cen = vec_length(a.x);
            diff = dist2cen - sizeL;
            if (diff > 0){
                // 1. mirror the position at the circular wall
                vec_mul221(wall_normal, 1 / dist2cen); // wall_normal = 1 * a.x / |a.x|
                hv = vec_mul(wall_normal, - 2 * diff); // hv = wall_normal * 2 * diff
                vec_add221(a.x, hv);
                // 2. mirror the velocity at the circular wall
                diff = vec_dot(wall_normal, a.v);
                hv = vec_mul(wall_normal, - 2 * diff);
                vec_add221(a.v, hv);
                // 3. update rest of agent properties
                a.u = vec_set_mag(a.v, 1);
                a.phi = atan2(a.u[1], a.u[0]);
            }
            break;
        case 6:
            // half-elastic circle boundary condition
            dist2cen = vec_length(a.x);
            diff = dist2cen - sizeL;
            if (diff > 0){
                // 1. mirror the position at the circular wall
                vec_mul221(wall_normal, 1 / dist2cen); // wall_normal = 1 * a.x / |a.x|
                hv = vec_mul(wall_normal, - 2 * diff); // hv = wall_normal * 2 * diff
                vec_add221(a.x, hv);
                // 2. set velocity component normal to wall to 0
                diff = vec_dot(wall_normal, a.v);
                hv = vec_mul(wall_normal, -diff);
                vec_add221(a.v, hv);
                vec_div221(a.v, 4);
                // 3. update rest of agent properties
                a.u = vec_set_mag(a.v, 1);
                a.phi = atan2(a.u[1], a.u[0]);
            }
            break;
    }
}

template
void Boundary(particle &a, double sizeL,  int BC);
