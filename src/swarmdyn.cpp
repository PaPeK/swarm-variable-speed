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
        }
        // if no run Sclu<MinCluster -> use one with largest cluster
        if (iter == maxiter - 1 && SysPara.Sclu < SysPara.MinCluster){
            agent = agentsave;
            SysPara.Sclu = bestSclu;
            // go back exactly where it stopped
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
    double rnv = 0.0;

#if LOCAL_METRIC==1
    int box=ptrSP->box;
#endif
    double dt = ptrSP->dt;
    int N = a.size();
    unsigned int j;
    int ii;

    // Reset simulation-step specific values to default
    if (N > 1){
        for (i=0; i<N; i++)
            a[i].NN.resize(0);
        InteractionVoronoiF2F(a, ptrSP);
    }

    // Update all agents
    for(i=0;i<N;i++)
    {
        // Generate noise
        rnp = ptrSP->noisep * gsl_ran_gaussian(r, 1.0);
        rnv = ptrSP->noisev * gsl_ran_gaussian(r, 1.0);
        MoveParticle(a[i], ptrSP, r, rnp, rnv);
    }
}


void Output(int s, std::vector<particle> &a, params &SP){
    std::vector<double> eddi = Dist2AlphaShape(a, &SP); // sets each agent.eddi
    std::vector<double> out;
    if (SP.out_mean){
        out = Out_swarm(a, SP);
        DataCreateSaveWrite(SP.dataOutMean, out, SP, "swarm");
    }
    if (SP.out_particle)
        WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
}


std::vector<double> Out_swarm(std::vector<particle> &a, params &SP){
    int N = a.size();
    double dist = 0;    // distance between predator and single prey
    std::vector<double> avg_x(2);   // average position vector
    // helper
    double hd = 0;
    std::vector<double> hv(2, 1);
    // to compute the aspect ratio:
    double max_IID = 0;
    std::vector<double> max_IID_vec(2, 1);
    // compute averages
    std::vector<int> allprey(N);
    std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
    std::vector<double> basic_avgs = basic_particle_averages(a, allprey);
    std::vector<double> avg_v {basic_avgs[0], basic_avgs[1]};
    std::vector<double> avg_u {basic_avgs[2], basic_avgs[3]};
    double avg_vsquare = basic_avgs[5];
    double avg_speed = basic_avgs[4];

    double IID = 0;     // Inter Individual Distance
    double NND = 0;     // Nearest Neighbor Distance
    double ND = 0;      // Neighbor Distance
    double nd = 0;
    double av_eddi = 0;
    double inv_NND2 = 0;
    std::vector<double> all_inv_NND2(N);
    for(int i=0; i<N; i++){
        vec_add221(avg_x, a[i].x);
        // NND:
        nd = 0;
        for (auto it=a[i].NN.begin(); it!=a[i].NN.end(); ++it){
            dist = CalcDist(a[i].x, a[*it].x, SP.BC, SP.sizeL);
            nd += dist;
            hd = fmin(hd, dist);
        }

        nd /= a[i].NN.size();
        ND += nd;

        // NND: (not in NN-loop because async-update do not has always NN)
        // ATTENTION: the computatoin of the NND uses both loops (also IID-loop)
        hd = SP.N * 10 * SP.rep_range;  // arbitrary large value
        for(int j=0; j<i-1; j++){
            dist = CalcDist(a[i].x, a[j].x, SP.BC, SP.sizeL);
            hd = fmin(hd, dist);
        }
        // IID:
        for(int j=i+1; j<N; j++){
            hv = CalcDistVec(a[i].x, a[j].x, SP.BC, SP.sizeL);
            dist = vec_length(hv);
            hd = fmin(hd, dist); // for NND
            IID += dist;
            if(dist > max_IID){ // maxIID and vector for aspect_ratio
                max_IID = dist;
                max_IID_vec = hv; 
            }
        }
        NND += hd;
        hd = 1 / (hd * hd);
        inv_NND2 += hd;
        all_inv_NND2[i] = hd;
        av_eddi += a[i].eddi;

    }

    vec_div221(avg_x, N);
    NND /= N;
    ND /= N;
    inv_NND2 /= N;
    IID /= (N-1) * N;
    av_eddi /= N;
    double avgvel = vec_length(avg_v);
    double pol_order = vec_length(avg_u);
    // double avgdirection = atan2(avg_v[1], avg_v[0]);

    //VARIANCE computation:
    double var_inv_NND2 = 0;
    for(int i=0; i<N; i++){
        hd = all_inv_NND2[i] - inv_NND2;
        var_inv_NND2 += hd * hd;
    }
    var_inv_NND2 /= N;

    // Calc normalized angular momentum (normalized by radial distance)
    // milling OP:
    double L_norm = 0;
    for(int i=0; i<N; i++){
        hv[0] = a[i].x[0] - avg_x[0];
        hv[1] = a[i].x[1] - avg_x[1];
        L_norm += (hv[0] * a[i].v[1] - hv[1] * a[i].v[0]) / (vec_length(hv));   // L/r=(\vec{r} x \vec{v})/r
    }
    L_norm = fabs(L_norm) / N;

    double Area_ConHull = AreaConvexHull(a, allprey);
    // computes elongation and aspect ratio:
    double elongation = get_elongation(a, avg_u, allprey);
    double aspect_ratio = get_elongation(a, max_IID_vec, allprey);
    double a1, a2, a_com_maxIID;
    a1 = (avg_v[0] * max_IID_vec[0] +
          avg_v[1] * max_IID_vec[1]) /
         (avgvel * vec_length(max_IID_vec));
    a2 = acos(-a1);
    a1 = acos(a1);
    a_com_maxIID = fmin(a1, a2); 
    // computes direction-correlation fctn.
    // double C_velFlucDir = corr_velFlucDirection_interact_all(a, avg_v);
    double varVelFluc = get_varOFVecFluctuation(a, avg_v, &particle::v, allprey);
    double C_velFluc = corr_velFluc_interact(a, avg_v);


    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(21);
    out_vec.push_back( N); 
    out_vec.push_back( pol_order); 
    out_vec.push_back( L_norm); 
    out_vec.push_back( avg_x[0]); 
    out_vec.push_back( avg_x[1]); 
    out_vec.push_back( avg_v[0]); 
    out_vec.push_back( avg_v[1]); 
    out_vec.push_back( avg_speed); 
    out_vec.push_back( avg_vsquare); 
    out_vec.push_back( elongation); 
    out_vec.push_back( aspect_ratio); 
    out_vec.push_back( a_com_maxIID); 
    out_vec.push_back( Area_ConHull); 
    out_vec.push_back( NND); 
    out_vec.push_back( IID); 
    out_vec.push_back( ND); 
    out_vec.push_back( varVelFluc);
    out_vec.push_back( C_velFluc); 
    out_vec.push_back( av_eddi);
    out_vec.push_back( inv_NND2);
    out_vec.push_back( var_inv_NND2);

    return out_vec;
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
            if ( SP.out_extend ){    // multiple Runs -> open and extend existing files
                h5dset = new H5::DataSet(H5out->openDataSet(n_dset.c_str()));
                h5read_extend_dim(h5dset, dim);
            }
            else{
                dim[0] = SP.total_outstep - SP.outstep; // is less/equal for particleD
                // dim[1] = a.size();  // could be problematic if agents already killed -> BUG
                dim[1] = Nori;   // avoids problem with killed prey (id might be larger than a.size)
                dim[2] = out.size();
                h5dset = new H5::DataSet(h5CreateDSet(H5out, dim,
                                                      n_dset.c_str(), "double"));
            }
        }
        // else: load dataSet and read dimension
        else{
            h5dset = new H5::DataSet(H5out->openDataSet(n_dset.c_str()));
            h5readDimension(h5dset, dim);
        }
        // now create-offset
        std::vector<hsize_t> offset(dim.size());  // offset of hyperslab  
        if (SP.out_extend)
            offset[0] = dim[0]-1;  // to not overwrite preceding run
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


template<class T>
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<T> &nodes){
    int N = nodes.size();
    //initialize
    unsigned int ii = 0;
    std::vector<double> avg_v(2), avg_u(2);
    double vsquare, avg_s, avg_vsquare;
    avg_s = avg_vsquare = 0;
    // compute output
    for(int i=0; i<N; i++){
        ii = nodes[i];
        vec_add221(avg_v, a[ii].v);
        vec_add221(avg_u, a[ii].u);
        vsquare = a[ii].v[0] * a[ii].v[0] + a[ii].v[1] * a[ii].v[1];
        avg_s += sqrt(vsquare);
        avg_vsquare += vsquare;
    }
    vec_div221(avg_v, N);
    vec_div221(avg_u, N);
    avg_s /= N;
    avg_vsquare /= N;
    // return output
    std::vector<double> out {avg_v[0], avg_v[1], avg_u[0], avg_u[1], avg_s, avg_vsquare};
    return out;
}
template
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<int> &nodes);
template
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<unsigned int> &nodes);


// Interesting: function takes member specifier
// function-call: get_var_of_vecNorm(a, mean, &particle::v) for velocity vector v
template<class T, class MembType>
double get_varOFVecFluctuation(std::vector<particle> &a,
                                 std::vector<double> &mean,
                                 MembType vec,
                                 std::vector<T> &nodes){
    int N = nodes.size();
    unsigned int ii;
    double var = 0;
    std::vector<double> v;
    for(int i=0; i<N; i++){
        ii = nodes[i];
        v = a[ii].*vec;
        vec_sub221(v, mean);
        // var += a[ii].*vec[0] * a[ii].*vec[0] + a[ii].*vec[1] * a[ii].*vec[1];
        var += v[0] * v[0] + v[1] * v[1];
    }
    var /= N;
    return var; 
}
template
double get_varOFVecFluctuation(std::vector<particle> &a,
                                 std::vector<double> &mean,
                                 std::vector<double> particle::* vec,
                                 std::vector<int> &nodes);
template
double get_varOFVecFluctuation(std::vector<particle> &a,
                                 std::vector<double> &mean,
                                 std::vector<double> particle::* vec,
                                 std::vector<unsigned int> &nodes);


// in comparison with corr_velFluc_interact, it evaluates all pair-correlations
// between interacting agents
// instead of only correlations between pairs belonging to specific sets
double corr_velFluc_interact(std::vector<particle> &a,
                                          std::vector<double> &mean_v){
    double corr_velFlucDir = 0;
    int N_interPairs = 0;   // # of interacting Fnn-Fnn2 pairs
    for(unsigned int i=0; i<a.size(); i++){
        std::vector<double> dv_i = vec_sub(a[i].v, mean_v);
        for(unsigned int j=0; j<a[i].NN.size(); j++){
            unsigned int jj = a[i].NN[j];
            std::vector<double> dv_j = vec_sub(a[jj].v, mean_v);
            corr_velFlucDir += vec_dot(dv_i, dv_j);
            N_interPairs++;
        }
    }
    if (N_interPairs > 0)
        corr_velFlucDir /= N_interPairs;
    return corr_velFlucDir;
}


void DataCreateSaveWrite(std::vector< std::vector<double> > &data,
                     std::vector<double> &out, params &SP,
                     std::string name, bool forceSave){
    unsigned int n = data.size();
    unsigned int m = out.size();
    // CREATE (initialize data-vector):
    if (n == 0){
        unsigned int left_outsteps = SP.total_outstep - SP.outstep;
        std::vector< std::vector<double> > vec2d(left_outsteps, std::vector<double>(m));
        data = vec2d;
        n = left_outsteps;
    }
    // SAVE current values to data-array
    unsigned int current = SP.outstep - (SP.total_outstep - n);
    data[current] = out;
    // WRITE TO FILE (last output-step)
    if ( (SP.outstep == SP.total_outstep - 1) || forceSave ){
        if (SP.out_h5){
            std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
            h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
        }
        else
            WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                          data, false);
    }
}
