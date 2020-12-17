/*  Settings
    defines functions to set parameters for
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Romanczuk

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
#include "settings.h"
void SetCoreParameters(params* SP)
{
    // Function setting core system parameters to prevent crashes
    SP->sizeL=100.;
    SP->N=100;
    SP->sim_time=100;
    SP->dt=0.05;
    SP->noisep=0.01;
    SP->sim_steps=2000;
    SP->trans_time=0;

    SP->speed0=1.0;

    SP->rep_range=1.0;
    SP->rep_strength=0.0;
    SP->rep_steepness=-0.2;

    SP->alg_strength=0.0;

    SP->output=1.0;
    SP->BC=-1;
    SP->IC=0;

    SP->step_output=1;
    SP->trans_time=0.0;

    SP->cludist = 5;   
    SP->Sclu = 0;
    SP->MinCluster = 1;
    SP->location = "./";
    SP->fileID = "xx";
    SP->out_h5 = 1;
    SP->outstep = 0;
}

void InitSystemParameters(params* SP)
{
    // Function initializing system parameters by updating the
    // parameter data after parsing the command line

    SP->noisep = sqrt(SP->dt * 2 * SP->Dphi);

    // Set output step
    if(SP->output < SP->dt)
        SP->output = SP->dt;
    SP->step_output=(int) (SP->output/SP->dt);
    // Set trans_time
    if(SP->trans_time >= SP->sim_time) // always at least 1 outputstep -> check if clustersize >= MinClusterSize
        SP->trans_time = SP->sim_time - SP->output;
    if(SP->trans_time < 0)
        SP->trans_time = 0;
    }

    // if BC not periodic -> set systemsize laarge
    if (SP->BC < 0)
        SP->sizeL =  SP->N * SP->N;

    // output parameters:
    SP->total_outstep = (int)((SP->sim_steps - 1) / SP->step_output) - 
                        (int)(SP->trans_time / SP->dt) / SP->step_output +
                        (0 == fmod(SP->trans_time / SP->dt / SP->step_output, 1));   // out starts at s<=trans_time/dt
    if (SP->total_outstep < 0)
        SP->total_outstep = 0;
}


void InitSystem(std::vector<particle> &a, params SP)
{
    unsigned int i;
    // Set arrays to default values
    for(i=0; i<a.size(); i++){
        a[i].x.assign(2, 0.);
        a[i].v.assign(2, 0.);
        a[i].u.assign(2, 0.);
        a[i].force.assign(2, 0.);
        a[i].force_rep.assign(2, 0.);
        a[i].force_alg.assign(2, 0.);
        a[i].vproj=SP.speed0;
        a[i].counter_rep=0;
        a[i].id = i;
        a[i].alg_strength = SP.alg_strength;
    }
}

void ResetSystem(std::vector<particle> &a, params *ptrSP, bool out, gsl_rng *r)
{
    // Function resetting all variables and setting the initial conditions

    int i, ii;
    double  theta=0.0;
    double  sizeL=ptrSP->sizeL;
    int     IC=ptrSP->IC;
    int     N=ptrSP->N;
    double  speed0=ptrSP->speed0;
    double  tmpspeed=0.0;
    double  range = ptrSP->rep_range;

    // set initialization distance between particles
    range = 1 * ptrSP->rep_range;

    // Different definitions of initial conditions
    if(IC==1) // Directionally ordered, Spatially disordered
    {
        theta=2*M_PI*gsl_rng_uniform(r);
        for(i=0;i<N;i++)
            {
                a[i].x[0]=sizeL*gsl_rng_uniform(r);
                a[i].x[1]=sizeL*gsl_rng_uniform(r);
                a[i].phi=theta;
            }
    }
    else if(IC==99) // Read from a file
    {
        LoadCoordinates(ptrSP, "initCoord", a, N, sizeL);
        for(i=0;i<N;i++)
            {
                tmpspeed=vec_length(a[i].v);
                if(tmpspeed>0.0)
                {
                    a[i].u[0]=a[i].v[0]/tmpspeed;
                    a[i].u[1]=a[i].v[1]/tmpspeed;
                }
                a[i].phi=atan2(a[i].v[1],a[i].v[0]);
                a[i].vproj=tmpspeed;
            }
    }
    else if(IC==0)
    {
        for(i=0;i<N;i++) // Directionally, Spatially disordered in a square
            {
            theta=2*M_PI*gsl_rng_uniform(r);
            a[i].x[0]=fmin(sizeL,sqrt(N)*range)*gsl_rng_uniform(r);
            a[i].x[1]=fmin(sizeL,sqrt(N)*range)*gsl_rng_uniform(r);
            a[i].phi=theta;
        }

    }
    else if(IC==2)
    {
        theta=2*M_PI*gsl_rng_uniform(r);
        for(i=0;i<N;i++) // Directionally ordered, Spatially disordered in a square
            {
            a[i].x[0]=fmin(sizeL,sqrt(N)*range)*gsl_rng_uniform(r);
            a[i].x[1]=fmin(sizeL,sqrt(N)*range)*gsl_rng_uniform(r);
            a[i].phi=theta;
        }

    }
    else if(IC==3) // Directionally ordered, Spatially ordered in circle but perturbed
    {
        theta=2*M_PI*gsl_rng_uniform(r);
        int circle = 1; // current circle
        int Ncir = 3;  // # of particles fitting on current circle (1st. circle less)
        int Ncap = Ncir;    // # of particle who does totally fit on all cicle including current one
        int ccircle = 0;    // # of particles on current circle
        double sigma_space = 0.8;
        // shuffle-order: to not always have a[0] in middle and a[-1] at periphery
        std::vector<int> allprey(N);
        std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
        unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::shuffle(allprey.begin(), allprey.end(), std::default_random_engine(seed));
        for(i=0; i<N; i++)
        {
            ii = allprey[i];
            if (i == Ncap){
                circle += 1;
                Ncir = int(2*M_PI*circle); // length of circumference
                Ncap += Ncir;
                ccircle = 0;
            }
            a[ii].x[0] = range*circle*cos(2*M_PI*ccircle/Ncir) 
                        + sigma_space*range*gsl_rng_uniform(r);
            a[ii].x[1] = range*circle*sin(2*M_PI*ccircle/Ncir) 
                        + sigma_space*range*gsl_rng_uniform(r);
            a[ii].phi=theta;
            a[ii].u[0] = cos(theta);
            a[ii].u[1] = sin(theta);
            a[ii].v[0] = a[ii].u[0];
            a[ii].v[1] = a[ii].u[1];
            ccircle += 1;
        }

    }
    else if(IC == 4) // Directionally milling-order, Spatially ordered in circle but perturbed
    {
        int circle = 3; // current circle
        int Ncir = int(2*M_PI*circle);  // # of particles fitting on current circle (1st. circle less)
        int Ncap = Ncir;    // # of particle which totally fit on all cicles including current one
        int ccircle = 0;    // # of particles on current circle
        double sigma_space = 1;
        double angle;
        // shuffle-order: to not always have a[0] in middle and a[-1] at periphery
        std::vector<int> allprey(N);
        std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
        std::random_shuffle(allprey.begin(), allprey.end());
        for(i=0; i<N; i++)
        {
            ii = allprey[i];
            if (i == Ncap){
                circle += 1;
                Ncir = int(2*M_PI*circle); // length of circumference
                Ncap += Ncir;
                ccircle = 0;
            }
            angle = 2*M_PI*ccircle/Ncir;
            a[ii].x[0] = range*circle*cos(angle) 
                        + sigma_space*range*gsl_rng_uniform(r);
            a[ii].x[1] = range*circle*sin(angle) 
                        + sigma_space*range*gsl_rng_uniform(r);
            angle += M_PI/2.;
            a[ii].phi = fmod(angle, 2*M_PI);
            ccircle += 1;
        }

    }
    else if(IC == 5){ // Directionally milling-order, Spatially ordered in circle but perturbed
        int circle = 3; // current circle
        int Ncir = int(2*M_PI*circle);  // # of particles fitting on current circle (1st. circle less)
        int Ncap = Ncir;    // # of particle which totally fit on all cicles including current one
        int ccircle = 0;    // # of particles on current circle
        double sigma_space = 1;
        double angle;
        // shuffle-order: to not always have a[0] in middle and a[-1] at periphery
        std::vector<int> allprey(N);
        std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
        std::random_shuffle(allprey.begin(), allprey.end());
        for(i=0; i<N; i++)
        {
            ii = allprey[i];
            if (i == Ncap){
                circle += 1;
                Ncir = int(2*M_PI*circle); // length of circumference
                Ncap += Ncir;
                ccircle = 0;
            }
            angle = 2*M_PI*ccircle/Ncir;
            a[ii].x[0] = range*circle*cos(angle) 
                        + sigma_space*range*gsl_rng_uniform(r);
            a[ii].x[1] = range*circle*sin(angle) 
                        + sigma_space*range*gsl_rng_uniform(r);
            theta = 2 * M_PI * gsl_rng_uniform(r);
            a[ii].phi = theta;
            ccircle += 1;
        }

    }
    else // e.g. IC==4
    {
        for(i=0;i<N;i++) // Directionally, Spatially dis-ordered
            {
            theta=2*M_PI*gsl_rng_uniform(r);
            a[i].x[0]=sizeL*gsl_rng_uniform(r);
            a[i].x[1]=sizeL*gsl_rng_uniform(r);
            a[i].phi=theta;
        }

    }
    if (IC != 99){  // not read from file -> all unit speed
        for(i=0;i<N;i++)
            {
            theta = a[i].phi;
            a[i].v[0] = cos(theta);
            a[i].v[1] = sin(theta);
            a[i].u[0] = cos(theta);
            a[i].u[1] = sin(theta);
            a[i].vproj = speed0;
        }
    }
    // Write initial coordinates of agents
    if (out){
        std::string fname = "init_coord.dat";
        WritePosVel(a, ptrSP, fname);
    }
    if( ptrSP->BC != -1 )
        for(i=0; i<N; i++)
            Boundary(a[i], ptrSP->sizeL, ptrSP->BC);
}
