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
    std::vector<particle> agent_dead(0);
    std::vector<particle> dummy_dead(0);
    InitSystem(agent, SysPara);
    InitRNG();
    ResetSystem(agent, &SysPara, false, r);
    LoadCoordinatesCPP(&SysPara, "init_posvel_" + SysPara.fileID, agent);
    pava_load_set_paras(agent, SysPara.location + "pava_in_" +
                        SysPara.fileID + ".in");
    double dt = SysPara.dt;

    unsigned int maxiter = 3;
    unsigned int iter = 0;
    std::vector<int> best_clu;
    std::vector<int> fnn;
    unsigned int bestSclu = 0;
    int sstart = 0;
    std::vector<predator> predsave;
    std::vector<particle> agentsave;
    std::vector<double> memosave;
    double t1 = clock(), t2 = 0.; // time variables for measuring comp. time
    std::cout<< std::endl;
    unsigned int Ndead_before = 0;
    bool pred_state_was0 = false; // necessary because break-condition only checked during output
    while (SysPara.Sclu < SysPara.MinCluster && iter < maxiter){
        // Perform numerical integrate
        t1= clock();
        for(s=sstart; s < SysPara.sim_steps; s++){
            // define some basic time-flags
            bool time_output = (s >= static_cast<int>(SysPara.trans_time/dt));
            bool time_1stoutput = (s == static_cast<int>(SysPara.trans_time/dt));
            // only simulate largest cluster if recording starts or predator appears:
            if (time_1stoutput && SysPara.cluster_attacked.size() == 0){
                SysPara.cluster_attacked = GetLargestCluster(agent, &SysPara);
                SysPara.Sclu = SysPara.cluster_attacked.size();
            }
            // EITHER at 1.st output OR at predator apprearence (not at both occasions)
            if (time_predAppears || time_1stoutput)
                if (SysPara.Npred > 0 && agent_dead.size() == 0) // ensures no re-split
                    split_notInCluster(agent, agent_dead, SysPara.cluster_attacked, preds);
            if (time_predAppears){
                // CREATE PREDATOR
                if ( preds.size() == 1 )
                    CreatePredator(agent, &SysPara, preds[0]);
                // copy Fs and P to dummies
                if (SysPara.out_dummy){
                    dummy = agent;
                    dummy_dead = agent_dead;
                    predsD = preds;
                }
            }
            Ndead_before = SysPara.Ndead;
            // Perform a single step
            // first split: output handles agents who are dead but
            //  NN of non-dead agents (also imporant for NN2)
            split_dead(agent, agent_dead, preds);
            Step(s, agent, &SysPara, preds);
            if (time_pred && SysPara.out_dummy){
                split_dead(dummy, dummy_dead, predsD);
                Step(s, dummy, &SysPara, predsD, true);
            }
            if (time_pred && SysPara.outstep > 10 && preds[0].state == 0)
                pred_state_was0 = true;

            // Data output
            if(s%SysPara.step_output==0 && time_output)
            {
                // RESTART-CONDITIONS unless maxiteration reached -> take the best run:
                if (SysPara.outstep == 0 && iter < maxiter - 1){
                    // if S(Cluster) < MinCluster
                    if (SysPara.Sclu < SysPara.MinCluster){
                        merge_dead(agent, agent_dead);
                        if (SysPara.Sclu > bestSclu){
                            memosave = SysPara.pred_memo_av;
                            agentsave = agent;
                            predsave = preds;
                            best_clu = SysPara.cluster_attacked;
                            bestSclu = SysPara.Sclu;
                        }
                        SysPara.Sclu = 0;
                        break;
                    }
                }
                // BREAK AND REPEAT conditions
                if ( time_pred && ( preds.size() == 1 ) ){
                    // if "undesired" Output happend -> unvalid run (Sclu=0) -> repeat
                    // TODO: break conditions only if preds.size() == 1
                    predator pred = preds[0];   // designed for preds.size() == 1
                    if (time_predAppears){
                        fnn = GetPredFrontPrey<int>(agent, &SysPara, &pred, pred.cluster);
                        if (fnn.size() == 0){
                            std::cout<< "no prey in FRONT!!!!" << std::endl;
                            ErrorOutput(agent, 1, &SysPara, &pred);
                            SysPara.Sclu = 0;
                            merge_dead(agent, agent_dead);
                            iter--; // iteration only counted for Sclu < MinCluster criterion
                            break;
                        }
                    }
                    // Break-condition dependent on predator-movement-scheme
                    bool breakAtReturn = false;
                    if (fmod(SysPara.pred_attack, 10) < 2 && !breakAtReturn)
                        breakAtReturn = (pred_state_was0 && SysPara.outstep > 10);
                        // if pred_com < 2 it moves always straight after it passed the com
                        //      -> therfore no turning and therefore stops if no prey in front
                    // run break procedure (write output and break)
                    if ((SysPara.stopAtKill && SysPara.Ndead >= 1) || breakAtReturn){
                        Output(s, agent, SysPara, preds, dummy, predsD, true);
                        break;
                    }
                }
                Output(s, agent, SysPara, preds, dummy, predsD);
                SysPara.outstep += 1;
                if (time_pred)
                    SysPara.outstep_pred += 1;
            }
        }
        iter++;
        // only reset the system it run was "unsuccesfull" and if it is not the final-run
        if (SysPara.Sclu < SysPara.MinCluster && iter < maxiter){
            InitSystem(agent, SysPara);
            ResetSystem(agent, &SysPara, false, r);
            SysPara.outstep = 0;
            SysPara.outstep_pred = 0;
            sstart = 0;
            LoadCoordinatesCPP(&SysPara, "init_posvel_"
                               + SysPara.fileID, agent);
            pava_load_set_paras(agent, SysPara.location + "pava_in_" +
                            SysPara.fileID + ".in");
        }
        // if no run Sclu<MinCluster -> use one with largest cluster
        if (iter == maxiter - 1 && SysPara.Sclu < SysPara.MinCluster){
            agent = agentsave;
            dummy = agentsave;
            preds = predsave;
            predsD = predsave;
            SysPara.pred_memo_av = memosave;
            SysPara.Sclu = bestSclu;
            SysPara.cluster_attacked = best_clu;
            if (SysPara.Npred > 0)
                split_notInCluster(agent, agent_dead, SysPara.cluster_attacked, preds);
            // go back exactly where it stopped
            std::vector<double> eddi = Dist2AlphaShape(agent, &SysPara);
            Output(s, agent, SysPara, preds, dummy, predsD);
            SysPara.outstep += 1;
            if ((s >= static_cast<int>(SysPara.pred_time/dt)))
                SysPara.outstep_pred += 1;
            sstart = static_cast<int>(SysPara.trans_time/dt) + 1;
        }
    }
    std::cout << "iters: "<< iter <<" Sclu: "<< SysPara.Sclu 
              << " MinCluster: "<< SysPara.MinCluster << "\n"; 
    // if minimum output generated -> assumes equilibration run
    // -> save final positions velocities
    merge_dead(agent, agent_dead);
    if (SysPara.outstep == 1)
        WritePosVel(agent, &SysPara, "final_posvel_" + SysPara.fileID, false);
    pava_output(agent, SysPara.cluster_attacked,
                SysPara.location + "pava_out_" + SysPara.fileID + ".dat");
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

void Step(int s, std::vector<particle> &a, params* ptrSP, std::vector<predator> &preds, bool dummy)
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

    // CREATE PREDATOR MEMORY //////////////////////////////
    // Start of memory: estimate largest cluster for hunt
    if (s >= static_cast<int>((ptrSP->pred_time - ptrSP->pred_memo_t)/dt) &&
        s-1 < static_cast<int>((ptrSP->pred_time - ptrSP->pred_memo_t)/dt)){
        ptrSP->cluster_attacked = GetLargestCluster(a, ptrSP);
        ptrSP->Sclu = ptrSP->cluster_attacked.size();
        for (j=0; j<preds.size(); j++)
            preds[j].cluster = ptrSP->cluster_attacked;
    }
    // during memory: average velocity of swarm
    // ATTENTION: "s < ptrSP->pred_time/dt" is important because "s-1 <ptrSP->pred_time/dt"
    //      will cause segmentation fault if cluster_attacked.size<N
    //          REASON: accessing particles via indices in once saved cluster_attacked
    //                  but indices are different at time1st_output -> split_dead function applied
    if (s >= static_cast<int>((ptrSP->pred_time - ptrSP->pred_memo_t)/dt) && 
        s < static_cast<int>(ptrSP->pred_time/dt)){
        double avvelx = 0;                 // Average x velocity of particles
        double avvely = 0;                 // Average y velocity of particles
        for (j = 0; j < ptrSP->cluster_attacked.size(); j++)
        {
            ii = ptrSP->cluster_attacked[j];
            avvelx += a[ii].v[0];
            avvely += a[ii].v[1];
        }
        avvelx /= ptrSP->cluster_attacked.size();
        avvely /= ptrSP->cluster_attacked.size();
        ptrSP->pred_memo_av[0] += avvelx;
        ptrSP->pred_memo_av[1] += avvely;
    }
    // Reset simulation-step specific values to default
    for (i=0; i<N; i++)
        a[i].NN.resize(0);
    for (i=0; i<preds.size(); i++){
        preds[i].NNset.clear();
        preds[i].NN.resize(0);
    }

    if (s < ptrSP->pred_time/dt)
        InteractionVoronoiF2F(a, ptrSP);
    else
        InteractionVoronoiF2FP(a, ptrSP, preds, dummy); // has build in pred->voronoinn computation

    // Update all agents AND dummy-agents (with same noise)
    for(i=0;i<N;i++)
    {
        // Generate noise
        rnp = ptrSP->noisep * gsl_ran_gaussian(r, 1.0);
        bool var_speed = false;
        MoveParticle(a[i], ptrSP, r, rnp);
    }
    // PREDATOR RELATED STUFF(P-move, Fitness, .... )
    if (s>=ptrSP->pred_time/dt){
        std::vector<int> allprey(N);
        std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
        std::vector<unsigned int> ui_vec;   // indices of F considered  by P(noticed by interaction and in front)
        // UPDATE P and initialize min, max
        for (i=0; i<N; i++)
            a[i].fit_decrease = 0;
        for (i=0; i<preds.size(); i++){
            // OVERWRITE pred->NN (already computed in Interaction BUT predator and agents did move!)
            ui_vec = GetPredVoronoiNN(a, ptrSP, &preds[i], allprey);
            preds[i].NN = GetPredFrontPrey<unsigned int, unsigned int>(a, ptrSP, &preds[i], ui_vec);
            Fitness(a, &preds[i], ptrSP, ptrSP->kill_range, r);
            MovePredator(ptrSP, &preds[i], a, s);
        }
    }
}


void Output(int s, std::vector<particle> &a, params &SP, std::vector<predator> &preds,
            std::vector<particle> &d, std::vector<predator> &predsD, bool forceSave){
    std::vector<double> out;

    if (s < SP.pred_time / SP.dt){
        WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
    }
    else{
        WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
        WriteParticles<predator>(preds, SP, "pred", SP.outstep_pred, SP.Npred);
        if (SP.out_dummy){
            WriteParticles<particle>(d, SP, "partD", SP.outstep_pred, SP.N);
            WriteParticles<predator>(predsD, SP, "predD", SP.outstep_pred, SP.Npred);
        }
    }
}


void ErrorOutput(std::vector<particle> &a, int err, params *ptrSP, 
             predator *pred){
    FILE *fp;
    unsigned int i;
    unsigned int N = a.size();
    fp=fopen((ptrSP->location + "E" + std::to_string(err) + "particle" + ptrSP->fileID + ".dat").c_str(), "a");
    for(i=0; i<N; i++){
        fprintf(fp,"%.4f\t%.4f\t%.4f\t%.4f\n", a[i].x[0], a[i].x[1], a[i].v[0], a[i].v[1]);
    }
    fprintf(fp,"\n\n");
    fclose(fp);
    
    // save predator values
    fp=fopen((ptrSP->location + "E" + std::to_string(err) + "predator" + ptrSP->fileID + ".dat").c_str(), "a");
    
    fprintf(fp,"%.4f\t%.4f\t%.4f\t%.4f\n", pred->x[0], pred->x[1], pred->v[0], pred->v[1]);
    fprintf(fp,"\n\n");
    fclose(fp);
    
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
template
void WriteParticles(std::vector<predator> &a, params &SP, 
                    std::string name, double outstep, int Nori);
