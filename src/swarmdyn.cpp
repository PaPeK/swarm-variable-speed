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

    std::vector<predator>  preds(SysPara.Npred);
    InitPredator(preds);
    // Create Dummies:
    std::vector<predator> predsD = preds;         // predator who hunts for dummy
    std::vector<particle> dummy = agent;    // dummy particles or prey (not seeing Predator)

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
            bool time_pred = (s >= static_cast<int>(SysPara.pred_time/dt));
            bool time_predAppears = (s == static_cast<int>(SysPara.pred_time/dt));
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
            // if (SysPara.Ndead > Ndead_before)
            //     WritePosVelDead(agent, SysPara, "pos_vel_dead_", pred);
            if (time_pred && SysPara.outstep > 10 && preds[0].state == 0)
                pred_state_was0 = true;

            // Data output
            if(s%SysPara.step_output==0 && time_output)
            {
                // UNCOMMENT if simulation per output step are interesting:
                // t2= clock();
                // double tdiff=(t2-t1)/CLOCKS_PER_SEC;
                // printf("s=%d; time / out. step = %.4f\n",s,tdiff);
                // t1=t2;
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
                        // if pred_com < 4 it moves always straight after it passed the com
                        //      -> therfore no turning and therefore stops if no prey in front
                    // run break procedure (write output and break)
                    if ((SysPara.stopAtKill && SysPara.Ndead >= 1) || breakAtReturn){
                        Output(s, agent, SysPara, preds, dummy, predsD, true);
                        break;
                    }
                }
                // if run is valid -> save position, velocity
                if (SysPara.outstep == 0){
                    std::vector<double> eddi = Dist2AlphaShape(agent, &SysPara);
                    WritePosVelPlus(agent, eddi, &SysPara,
                                    "pos_vel_eddi_" + SysPara.fileID, false);
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
            // if (SysPara.out_mean || SysPara.out_particle)     // for evolution-simulation no output needed
            std::vector<double> eddi = Dist2AlphaShape(agent, &SysPara);
            WritePosVelPlus(agent, eddi, &SysPara,
                            "pos_vel_eddi_" + SysPara.fileID, false);
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

    if (SP.outstep == 0) // saves initial x0, x1, v0, v1, eddi
        Out_start_Save(a, "start", SP);
    if (s < SP.pred_time / SP.dt){
        if (SP.out_mean){
            out = Out_swarm(a, SP);
            DataCreateSaveWrite(SP.dataOutMean, out, SP,
                                "swarm", forceSave);
        }
        if (SP.out_particle)
            WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
    }
    else{
        if (SP.outstep_pred == 0) // Attention: only for Npred=1
            Out_startPred_Save(preds, "start_pred", SP);
        if (SP.out_mean){
            out = Out_swarm(d, SP);
            DataCreateSaveWrite(SP.dataOutSwarm, out, SP,
                                "swarm", forceSave);
            out = Out_swarm_pred(a, preds[0], SP);
            DataCreateSaveWrite(SP.dataOutSwarmPred, out, SP, 
                                "swarm_pred", forceSave);
            if (SP.out_dummy){
                out = Out_swarm_pred(d, predsD[0], SP);
                DataCreateSaveWrite(SP.dataOutSwarmPredD, out, SP, 
                                    "swarm_predD", forceSave);
            }
        }
        if (SP.out_particle){
            WriteParticles<particle>(a, SP, "part", SP.outstep, SP.N);
            WriteParticles<predator>(preds, SP, "pred", SP.outstep_pred, SP.Npred);
            if (SP.out_dummy){
                WriteParticles<particle>(d, SP, "partD", SP.outstep_pred, SP.N);
                WriteParticles<predator>(predsD, SP, "predD", SP.outstep_pred, SP.Npred);
            }
        }
    }
    if ((SP.outstep == SP.total_outstep - 1) || forceSave){
        Out_end_Save(a, "end", SP);
        if (preds.size() >= 1)
            Out_End_PosVel(a, preds[0], "end_PosVel", SP);
        if (SP.out_dummy)
            Out_end_Save(d, "endD", SP);
    }
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

    }
    vec_div221(avg_x, N);
    NND /= N;
    ND /= N;
    IID /= (N-1) * N;
    double avgvel = vec_length(avg_v);
    double pol_order = vec_length(avg_u);
    // double avgdirection = atan2(avg_v[1], avg_v[0]);

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
    double C_velFlucDir = corr_velFlucDirection_interact_all(a, avg_v);
    double C_velFluc = corr_velFluc_interact_all(a, avg_v);

    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(17);
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
    out_vec.push_back( C_velFlucDir); 
    out_vec.push_back( C_velFluc); 

    return out_vec;
}


std::vector<double> Out_swarm_pred(std::vector<particle> &a,
                                   predator &pred, params &SP){
    std::vector<double> basic_avgs;
    std::vector<int> allprey(a.size());
    std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
    // CLU  ############################################################
    int N_clu = SP.cluster_attacked.size();
    // averaging
    basic_avgs = basic_particle_averages(a, allprey);
    std::vector<double> clu_avg_v {basic_avgs[0], basic_avgs[1]};
    double clu_avg_s = basic_avgs[4];
    double clu_avg_vsquare = basic_avgs[5];
    double clu_dpi = get_avg_pred_dist(a, allprey, pred, SP);
    double fit_sum = 0;
    double fit_count = 0;
    unsigned int ii;
    for(int i=0; i<a.size(); i++){
        if (a[i].fitness < 0){
            fit_sum += a[i].fitness;
            fit_count += 1;
        }
    }
    // FKN  ############################################################
    std::vector<int> fkn = GetPredFrontPrey<int>(a, &SP, &pred, allprey);
    fkn = GetPredKnn(a, &SP, &pred, 10, fkn);
    int N_fkn = fkn.size();
    // averaging
    basic_avgs = basic_particle_averages(a, fkn);
    double fkn_avg_s = basic_avgs[4];
    double fkn_avg_vsquare = basic_avgs[5];
    double fkn_dpi = get_avg_pred_dist(a, fkn, pred, SP);
    // FNN  ############################################################
    std::vector<unsigned int> Fnn;     // contains indices of prey fronta l-NN of pred
    std::copy(pred.NNset.begin(), pred.NNset.end(), std::back_inserter(Fnn));
    Fnn = GetPredFrontPrey<unsigned int, unsigned int>(a, &SP, &pred, Fnn);
    int N_fnn = Fnn.size();
    // averaging
    double fnn_deg = get_avg_deg(a, Fnn);
    basic_avgs = basic_particle_averages(a, Fnn);
    double fnn_avg_s = basic_avgs[4];
    double fnn_avg_vsquare = basic_avgs[5];
    // FNN2 ############################################################
    find_NN2set(a, &pred);
    std::vector<unsigned int> Fnn2;
    std::copy(pred.NN2set.begin(), pred.NN2set.end(), std::back_inserter(Fnn2));
    Fnn2 = GetPredFrontPrey<unsigned int, unsigned int>(a, &SP, &pred, Fnn2);
    int N_fnn2 = Fnn2.size();
    // averaging
    basic_avgs = basic_particle_averages(a, Fnn2);
    double fnn2_avg_s = basic_avgs[4];
    double fnn2_avg_vsquare = basic_avgs[5];
    double fnn2_dpi = get_avg_pred_dist(a, Fnn2, pred, SP);
    std::vector<double> CFnnFnn2_velFlucDir2 = corr_velFluc_interact(a, clu_avg_v, Fnn, Fnn2);
    double CFnnFnn2_velFluc = CFnnFnn2_velFlucDir2[0];
    double CFnnFnn2_velFlucDir = CFnnFnn2_velFlucDir2[1];
    std::vector<double> corrFnnFnn2_velFluc_all = corr_velFlucTS_interact(a, clu_avg_v, Fnn, Fnn2);
    double CFnnFnn2_velFlucTS = corrFnnFnn2_velFluc_all[0];
    double CFnnFnn2_velFlucTSNormed = corrFnnFnn2_velFluc_all[1];
    double CFnnFnn2_speFlucTS = corrFnnFnn2_velFluc_all[2];
    double CFnnFnn2_speFlucTSNormed = corrFnnFnn2_velFluc_all[3];
    // averaging-raw
    double fnn2_qamp = 0; // (links to fnn nodes + links to fnn2 nodes)/degree
    double fnn2_qsig = 0; // (links to fnn nodes)/degree
    unsigned int jj;
    for(int i=0; i<N_fnn2; i++)
    {
        ii = static_cast<int>(Fnn2[i]);
        double k_fnn = 0;     // # of links to fnn
        double k_fnn2 = 0;    // # of links to fnn2
        for(unsigned int j=0; j<a[ii].NN.size(); j++){
            jj = a[ii].NN[j];
            if(Fnn2.end() != std::find(Fnn2.begin(), Fnn2.end(), jj))
                k_fnn2++;
            else if (Fnn.end() != std::find(Fnn.begin(), Fnn.end(), jj))
                k_fnn++;
        }
        fnn2_qamp += (k_fnn + k_fnn2) / a[ii].NN.size();
        fnn2_qsig += k_fnn / a[ii].NN.size();
    }
    fnn2_qamp /= N_fnn2;
    fnn2_qsig /= N_fnn2;
    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(28);
    out_vec.push_back( N_clu);
    out_vec.push_back( clu_avg_s);
    out_vec.push_back( clu_avg_vsquare);
    out_vec.push_back( clu_dpi);
    out_vec.push_back( fit_sum);
    out_vec.push_back( fit_count);
    out_vec.push_back( N_fkn);
    out_vec.push_back( fkn_avg_s);
    out_vec.push_back( fkn_avg_vsquare);
    out_vec.push_back( fkn_dpi);
    out_vec.push_back( N_fnn);
    out_vec.push_back( fnn_avg_s);
    out_vec.push_back( fnn_avg_vsquare);
    out_vec.push_back( fnn_deg);
    out_vec.push_back( N_fnn2);
    out_vec.push_back( fnn2_avg_s);
    out_vec.push_back( fnn2_avg_vsquare);
    out_vec.push_back( fnn2_dpi);
    out_vec.push_back( fnn2_qamp);
    out_vec.push_back( fnn2_qsig);
    out_vec.push_back( CFnnFnn2_velFluc);
    out_vec.push_back( CFnnFnn2_velFlucDir);
    out_vec.push_back( CFnnFnn2_velFlucTS);
    out_vec.push_back( CFnnFnn2_velFlucTSNormed);
    out_vec.push_back( CFnnFnn2_speFlucTS);
    out_vec.push_back( CFnnFnn2_speFlucTSNormed);
    out_vec.push_back( pred.kills);
    return out_vec;
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


void Out_end_Save(std::vector<particle> &a, std::string name, params &SP){
    // create data:
    unsigned int n = SP.N;
    unsigned int m = 2;      // fit, (not_dead && in_clu)
    std::vector< std::vector<double> > data(n, std::vector<double>(m));
    for(int i=0; i < a.size(); i++){
        data[a[i].id][0] = a[i].fitness;
        data[a[i].id][1] = 1;
    }
    // WRITE TO FILE
    if (SP.out_h5){
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
    }
    else
        WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                      data, false);
}


void Out_End_PosVel(std::vector<particle> &a, predator &p, std::string name, params &SP){
    // purpose: get position and velocity data at the 
    //          end of the simulation to see if the predator leaves 
    //          many in front of him untouched
    // OUTPUT: 
    //      pv.shape(N+1, 4)
    //      position and velocity of all agents and the predator
    // create data:
    unsigned int n = SP.N + 1; // agents + 1 pred
    unsigned int m = 4;      // x0, x1, v0, v1
    unsigned int ii = 0;
    std::vector< std::vector<double> > data(n, std::vector<double>(m));
    for(int i=0; i < a.size(); i++){
        ii = a[i].id;
        data[ii][0] = a[i].x[0];
        data[ii][1] = a[i].x[1];
        data[ii][2] = a[i].v[0];
        data[ii][3] = a[i].v[1];
    }
    data[n-1][0] = p.x[0];
    data[n-1][1] = p.x[1];
    data[n-1][2] = p.v[0];
    data[n-1][3] = p.v[1];
    // WRITE TO FILE
    if (SP.out_h5){
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
    }
    else
        WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                      data, false);
}


void Out_start_Save(std::vector<particle> &a, std::string name, params &SP){
    // create data:
    unsigned int n = SP.N;
    unsigned int m = 5;      // x0, x1, v0, v1, eddi
    std::vector< std::vector<double> > data(n, std::vector<double>(m));
    std::vector<double> eddi = Dist2AlphaShape(a, &SP);
    unsigned int ii = 0;
    for(int i=0; i < a.size(); i++){
        ii = a[i].id;
        data[ii][0] = a[i].x[0];
        data[ii][1] = a[i].x[1];
        data[ii][2] = a[i].v[0];
        data[ii][3] = a[i].v[1];
        data[ii][4] = eddi[i];
    }
    // WRITE TO FILE
    if (SP.out_h5){
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
    }
    else
        WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                      data, false);
}


void Out_startPred_Save(std::vector<predator> &p, std::string name, params &SP){
    // create data:
    unsigned int n = p.size();
    std::vector<double> out = p[0].out();
    unsigned int m = out.size();
    std::vector< std::vector<double> > data(n, std::vector<double>(m));
    for(int i=0; i < p.size(); i++)
        data[i] = p[i].out();
    // WRITE TO FILE
    if (SP.out_h5){
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
    }
    else
        WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                    data, false);
}


/*Write_out: takes a std::vector and writes its content to 
* a hdf5-dataset (h5dset) or to a file (name) depending
* on the output mode (SP.out_h5)
*/
void Write_out(std::vector<double> &out, params &SP,
           std::vector<hsize_t> & vec_dim,
           H5::DataSet *h5dset, std::string name){
    if (SP.out_h5){
        // create offset of data 
        std::vector<hsize_t> offset(vec_dim.size());  // offset of hyperslab  
        if (SP.out_extend)
            offset[0] = vec_dim[0] - 1;  // to not overwrite preceding run
        offset[offset.size() - 2] = SP.outstep;      // corresponds to time offset
        h5WriteDouble(h5dset, out, offset);
    }
    else
        WriteVector(SP.location + name + "_" + SP.fileID + ".dat", out, false);
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

template<class T>
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<T> &nodes, predator &pred,
                         params &SP){
    T N = nodes.size();
    unsigned int ii;
    double dist = 0;
    double dpi;
    std::vector<double> r_pi(2);
    for(T i=0; i<N; i++){
        ii = nodes[i];
        r_pi = CalcDistVec(pred.x, a[ii].x, SP.BC, SP.sizeL);
        dpi = vec_length(r_pi);
        dist += dpi;
    }
    dist /= N;
    return dist; 
}
template
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<int> &nodes, predator &pred,
                         params &SP);
template
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<unsigned int> &nodes, predator &pred,
                         params &SP);

template<class T>
double get_avg_deg(std::vector<particle> &a, std::vector<T> &nodes){
    int N = nodes.size();
    unsigned int ii;
    double deg = 0;
    for(int i=0; i<N; i++){
        ii = nodes[i];
        deg += a[ii].NN.size();
    }
    deg /= N;
    return deg; 
}
template
double get_avg_deg(std::vector<particle> &a, std::vector<int> &nodes);
template
double get_avg_deg(std::vector<particle> &a, std::vector<unsigned int> &nodes);

/* -correlates velocity fluctuation in between 2 different sets of 
 *  agents (e.g.: Fnn=Transmitter, Fnn2=Receiver)
 * -ATTENTION: sets must be interacting (otherwise 0 correlation)
 * -since there is no clear definition of velocity correlations for vectors with an length != 1
 *  the correlation of the magnitude of the vecotors is used and additionally weighted with the
 *  cosine of the angle between the vectors
 * -thus the magnitude is kind of normalized by 
 *  switching direction of vector if it is supposed to be negative
 */
std::vector<double> corr_velFlucTS_interact(std::vector<particle> &a,
                                          std::vector<double> &mean_v,
                                          std::vector<unsigned int> Fnn,
                                          std::vector<unsigned int> Fnn2){
    // helper
    unsigned int ii, jj;
    std::vector<double> vec_h;
    double var_h, var_h2;
    // transform into a "timeseries" of agents
    std::vector<int> fnns, fnn2s;
    std::vector<double> fnns_dspe, fnn2s_dspe;
    fnns.reserve(2 * Fnn2.size());
    fnn2s.reserve(2 * Fnn2.size());
    fnns_dspe.reserve(2 * Fnn2.size());
    fnn2s_dspe.reserve(2 * Fnn2.size());
    double av_mag_fnn = 0;
    double av_mag_fnn2 = 0;
    int N_interPairs = 0;   // # of interacting Fnn-Fnn2 pairs
    for(unsigned int i=0; i<Fnn.size(); i++){
        ii = Fnn[i];
        for(unsigned int j=0; j<a[ii].NN.size(); j++){
            jj = a[ii].NN[j];
            vec_h = vec_sub(a[ii].v, mean_v);
            var_h = vec_length(vec_h);
            if (Fnn2.end() != std::find(Fnn2.begin(), Fnn2.end(), jj)){
                // fnn-partner
                fnns.push_back(ii);
                av_mag_fnn += var_h;
                fnns_dspe.push_back(var_h);
                // fnn2-partner
                fnn2s.push_back(jj);
                vec_h = vec_sub(a[jj].v, mean_v);
                var_h2 = vec_length(vec_h);
                av_mag_fnn2 += var_h2;
                fnn2s_dspe.push_back(var_h2);
                N_interPairs++;
            }
        }
    }
    av_mag_fnn /= N_interPairs;
    av_mag_fnn2 /= N_interPairs;
    double std_mag_fnn = 0;
    double std_mag_fnn2 = 0;
    double corr_speFluc = 0;
    double corr_velFluc = 0;
    for(unsigned int i=0; i<fnns.size(); i++){
        ii = fnns[i];
        jj = fnn2s[i];
        var_h = fnns_dspe[i] - av_mag_fnn;
        var_h2 = fnn2s_dspe[i] - av_mag_fnn2;
        std_mag_fnn += var_h * var_h;
        std_mag_fnn2 += var_h2 * var_h2;
        // speed-correlation:
        corr_speFluc += var_h * var_h2;
        // velocity-correlation:
        std::vector<double> dv_i = vec_sub(a[ii].v, mean_v);
        std::vector<double> dv_j = vec_sub(a[jj].v, mean_v);
        dv_i = vec_set_mag(dv_i, var_h);
        dv_j = vec_set_mag(dv_j, var_h2);
        corr_velFluc += vec_dot(dv_i, dv_j);
    }
    std_mag_fnn /= N_interPairs;
    std_mag_fnn2 /= N_interPairs;
    std_mag_fnn = sqrt(std_mag_fnn);
    std_mag_fnn2 = sqrt(std_mag_fnn2);
    corr_velFluc /= N_interPairs;
    corr_speFluc /= N_interPairs;
    double corr_velFlucNormed = corr_velFluc / (std_mag_fnn * std_mag_fnn2);
    double corr_speFlucNormed = corr_speFluc / (std_mag_fnn * std_mag_fnn2);
    std::vector<double> output {corr_velFluc, corr_velFlucNormed,
                                corr_speFluc, corr_speFlucNormed};
    return output;
}


std::vector<double> corr_velFluc_interact(std::vector<particle> &a,
                                                   std::vector<double> &mean_v,
                                                   std::vector<unsigned int> Fnn,
                                                   std::vector<unsigned int> Fnn2){
    // std::vector<unsigned int> Fnn;
    // std::copy(pred.NNset.begin(), pred.NNset.end(), std::back_inserter(Fnn));
    // Fnn = GetPredFrontPrey<unsigned int, unsigned int>(a, ptrSP, &pred, Fnn);
    double corr_velFlucDir = 0;
    double corr_velFluc = 0;
    int N_interPairs = 0;   // # of interacting Fnn-Fnn2 pairs
    for(unsigned int i=0; i<Fnn.size(); i++){
        unsigned int ii = Fnn[i];
        std::vector<double> dv_i = vec_sub(a[ii].v, mean_v);
        for(unsigned int j=0; j<a[ii].NN.size(); j++){
            unsigned int jj = a[ii].NN[j];
            if (Fnn2.end() != std::find(Fnn2.begin(), Fnn2.end(), jj)){
                std::vector<double> dv_j = vec_sub(a[jj].v, mean_v);
                double correlation = vec_dot(dv_i, dv_j);
                corr_velFluc += correlation;
                corr_velFlucDir += correlation / ( vec_length(dv_i) * vec_length(dv_j) );
                N_interPairs++;
            }
        }
    }
    corr_velFluc /= N_interPairs;
    corr_velFlucDir /= N_interPairs;
    std::vector<double> output {corr_velFluc, corr_velFlucDir};
    return output;
}


// in comparison with corr_velFluc_interact, it evaluates all pair-correlations
// between interacting agents
// instead of only correlations between pairs belonging to specific sets
double corr_velFlucDirection_interact_all(std::vector<particle> &a,
                                          std::vector<double> &mean_v){
    double corr_velFlucDir = 0;
    int N_interPairs = 0;   // # of interacting Fnn-Fnn2 pairs
    for(unsigned int i=0; i<a.size(); i++){
        std::vector<double> dv_i = vec_sub(a[i].v, mean_v);
        for(unsigned int j=0; j<a[i].NN.size(); j++){
            unsigned int jj = a[i].NN[j];
            std::vector<double> dv_j = vec_sub(a[jj].v, mean_v);
            corr_velFlucDir += vec_dot(dv_i, dv_j) / ( vec_length(dv_i) * vec_length(dv_j) );
            N_interPairs++;
        }
    }
    if (N_interPairs > 0)
        corr_velFlucDir /= N_interPairs;
    return corr_velFlucDir;
}


// in comparison with corr_velFluc_interact, it evaluates all pair-correlations
// between interacting agents
// instead of only correlations between pairs belonging to specific sets
double corr_velFluc_interact_all(std::vector<particle> &a,
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
        // if outstep == 0: create OR extend-dataset + read-dimension
        //      extend: easy... just extend 0 dimension
        //      no-extend: creat with dim(time=total_outstep-SP.outstep, N, out.size())
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
        // h5dset->close();
        // H5out->close();
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
