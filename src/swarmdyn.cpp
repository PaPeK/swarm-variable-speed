/*  SwarmDynamics 
    Stochastic differential equation model of agents interacting via attraction/repulsion/alignment.
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
#include "swarmdyn.h"

// GLOBAL STRUCT DECLARATION


// GLOBAL VARIABLE DECLARATION
long seed;          // seed rng

// gsl random number generator
gsl_rng *r; // global generator
const gsl_rng_type * T;

//Runs with ./swarmdyn
int main(int argc, char **argv){

    int s;              // variable for current step number eg: currenttime/dt

    params SysPara;   // system

    // Init and output the system paramters
    SetCoreParameters(&SysPara);
    ParseParameters(argc, argv, &SysPara);
    InitSystemParameters(&SysPara);
    OutputParameters(SysPara);

    // initialize agents and set initial conditions
    std::vector<particle> agent(SysPara.N);    // particles or prey
    InitSystem(agent, SysPara);
    InitRNG();
    ResetSystem(agent, &SysPara, false, r);
    LoadCoordinatesCPP(&SysPara, "init_posvel_" + SysPara.fileID, agent);
    pava_load_set_paras(agent, SysPara.location + "pava_in_" +
                        SysPara.fileID + ".in");
    double dt = SysPara.dt;

    unsigned int maxiter = 3;
    unsigned int iter = 0;
    std::vector<int> fnn;
    unsigned int bestSclu = 0;
    int sstart = 0;
    std::vector<particle> agentsave;
    double t1 = clock(), t2 = 0.; // time variables for measuring comp. time
    std::cout<< std::endl;
    while (SysPara.Sclu < SysPara.MinCluster && iter < maxiter){
        // Perform numerical integrate
        t1= clock();
        for(s=sstart; s < SysPara.sim_steps; s++){
            // define some basic time-flags
            bool time_output = (s >= static_cast<int>(SysPara.trans_time/dt));
            bool time_1stoutput = (s == static_cast<int>(SysPara.trans_time/dt));
            // compute the largest cluster:
            if (time_1stoutput){
                std::vector<int> largestCluster = GetLargestCluster(agent, &SysPara);
                SysPara.Sclu = largestCluster.size();
            }
            // Perform a single step
            Step(s, agent, &SysPara);
            // Data output
            if((s % SysPara.step_output == 0) and time_output)
            {
                // RESTART-CONDITIONS unless maxiteration reached -> take the best run:
                if (SysPara.outstep == 0 && iter < maxiter - 1){
                    // if S(Cluster) < MinCluster
                    if (SysPara.Sclu < SysPara.MinCluster){
                        if (SysPara.Sclu > bestSclu){
                            agentsave = agent;
                            bestSclu = SysPara.Sclu;
                        }
                        SysPara.Sclu = 0;
                        break;
                    }
                }
                Output(s, agent, SysPara);
                SysPara.outstep += 1;
            }
        }
        iter++;
        // only reset the system it run was "unsuccesfull" and if it is not the final-run
        if (SysPara.Sclu < SysPara.MinCluster && iter < maxiter){
            InitSystem(agent, SysPara);
            ResetSystem(agent, &SysPara, false, r);
            SysPara.outstep = 0;
            sstart = 0;
            LoadCoordinatesCPP(&SysPara, "init_posvel_"
                               + SysPara.fileID, agent);
            pava_load_set_paras(agent, SysPara.location + "pava_in_" +
                            SysPara.fileID + ".in");
        }
        // if no run Sclu<MinCluster -> use one with largest cluster
        if (iter == maxiter - 1 && SysPara.Sclu < SysPara.MinCluster){
            agent = agentsave;
            SysPara.Sclu = bestSclu;
            // go back exactly where it stopped
            std::vector<double> eddi = Dist2AlphaShape(agent, &SysPara);
            Output(s, agent, SysPara);
            SysPara.outstep += 1;
            sstart = static_cast<int>(SysPara.trans_time/dt) + 1;
        }
    }
    std::cout << "iters: "<< iter <<" Sclu: "<< SysPara.Sclu 
              << " MinCluster: "<< SysPara.MinCluster << "\n"; 
    // if minimum output generated -> assumes equilibration run
    // -> save final positions velocities
    if (SysPara.outstep == 1)
        WritePosVel(agent, &SysPara, "final_posvel_" + SysPara.fileID, false);
    return 0;
}


long unsigned int getseed(int const K)
{

    typedef std::chrono::high_resolution_clock hiclock;

    auto gett= [](std::chrono::time_point<hiclock> t0)
    {
        auto tn = hiclock::now();
        return static_cast<long unsigned int>(std::chrono::duration_cast<std::chrono::microseconds>(tn-t0).count());
    };

    long unsigned int diffs[10];
    diffs[0] = gett(hiclock::now());
    for(int i=1; i!=10; i++)
    {
        auto last = hiclock::now();
        for(int k=K; k!=0; k--)
        {
            diffs[i]= gett(last);
        }
    }

    return *std::max_element(&diffs[1],&diffs[9]);
}


void InitRNG(){
    // Initialize random number generator
    // time_t  t1;                     // Get system time for random number seed
    // time(&t1);
    // seed=t1;
    // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    // seed=std::chrono::time_cast<long int>(t1);
    seed = static_cast<unsigned long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
#if PRINT_PARAMS
    printf("Random number seed: %d\n",static_cast<int>(seed));
#endif
    gsl_rng_env_setup();
        T=gsl_rng_default;
        r=gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
}

void Step(int s, std::vector<particle> &a, params* ptrSP)
{
    // function for performing a single (Euler) integration step
    int i = 0;
    double rnp = 0.0;

#if LOCAL_METRIC==1
    int box=ptrSP->box;
#endif
    double dt = ptrSP->dt;
    int N = a.size();
    unsigned int j;
    int ii;

    // Reset simulation-step specific values to default
    for (i=0; i<N; i++)
        a[i].NN.resize(0);
    InteractionVoronoiF2F(a, ptrSP);

    // Update all agents
    for(i=0;i<N;i++)
    {
        // Generate noise
        rnp = ptrSP->noisep * gsl_ran_gaussian(r, 1.0);
        bool var_speed = false;
        MoveParticle(a[i], ptrSP, r, rnp);
    }
}


void Output(int s, std::vector<particle> &a, params &SP){
    // TODO: add more fctns for output
    WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
}


template<class part>
void WriteParticles(std::vector<part> &a, params &SP, 
                    std::string name, double outstep, int Nori){
    std::vector<double> out = a[0].out();
    if (SP.out_h5){
        H5::DataSet *h5dset;
        // load/create h5-file
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        H5::H5File *H5out;
        if ( exists(f_h5out) )
            H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_RDWR);
        else
            H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_TRUNC);
        std::vector<hsize_t> dim(3, 0);
        std::string n_dset;
        if (SP.fileID == "xx")
            n_dset = "/" + name;
        else
            n_dset = "/" + SP.fileID + "/" + name;
        if ( outstep == 0){
            dim[0] = SP.total_outstep - SP.outstep; // is less/equal for particleD
            // dim[1] = a.size();  // could be problematic if agents already killed -> BUG
            dim[1] = Nori;   // avoids problem with killed prey (id might be larger than a.size)
            dim[2] = out.size();
            h5dset = new H5::DataSet(h5CreateDSet(H5out, dim,
                                                  n_dset.c_str(), "double"));
        }
        // else: load dataSet and read dimension
        else{
            h5dset = new H5::DataSet(H5out->openDataSet(n_dset.c_str()));
            h5readDimension(h5dset, dim);
        }
        // now create-offset
        std::vector<hsize_t> offset(dim.size());  // offset of hyperslab  
        offset[dim.size()-3] = outstep;  // time-offset
        // write data to offset
        for(int i=0; i<a.size(); i++){
            out = a[i].out();
            offset[offset.size()-2] = a[i].id;  // corresponds to particle offset
            h5WriteDouble(h5dset, out, offset);
        }
        delete h5dset;
        delete H5out;
    }
    else {
        std::ofstream outFile((SP.location + name + "_" + SP.fileID
                               + ".dat").c_str(), std::ios::app);
        std::vector<double> out_default(out.size(), 0);
        unsigned int id = 0;
        for(int i=0; i<a.size(); i++){
            out = a[i].out();
            while (id < a[i].id){
                for (int j=0; j<out_default.size(); j++)
                    outFile << out_default[j] << " ";
                outFile << std::endl;
                id++;
            }
            for (int j=0; j<out.size(); j++)
                outFile << out[j] << " ";
            outFile << std::endl;
            id++;
        }
        while (id < Nori){
            for (int j=0; j<out_default.size(); j++)
                outFile << out_default[j] << " ";
            outFile << std::endl;
            id++;
        }
    }
}

template
void WriteParticles(std::vector<particle> &a, params &SP, 
                    std::string name, double outstep, int Nori);
