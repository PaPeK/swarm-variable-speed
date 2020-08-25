/*  AgentsDynamics
    defines movement rules/differential equations of agents in
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser
*/ 

#include "agents_dynamics.h"

std::vector<double> Prob_hit(std::vector<double> predNNdists, double kill_range){
    double dist;
    unsigned int j;
    std::vector<double> prob_hit(predNNdists.size(), 0);
    for (j=0; j<predNNdists.size(); j++){
        dist = predNNdists[j];
        if (dist < kill_range)
            prob_hit[j] = (kill_range - dist) / kill_range;
    }
    return prob_hit;
}


void Fitness(std::vector<particle> &a, predator *pred, params * ptrSP,
                      double kill_range, gsl_rng *r){
    if (pred->NN.size() == 0)
        return;
    double lucky;      // random number: lucky < danger -> killed
    double prob_catch;     // represent how likely prey gets killed in timestep dt (= -fitness_rate * dt)
    double rate_catch = ptrSP->kill_rate; // how many fish can it catch in 1 time unit
    unsigned int j, ii;
    double dist;
    std::vector<double> prob_hit;
    std::vector<double> prob_selected(pred->NN.size(), 0);
    // compute distances (multiple times needed)
    std::vector<double> predNNdists(pred->NN.size());
    for (j=0; j<pred->NN.size(); j++){
        ii = pred->NN[j];
        dist = CalcDist(pred->x, a[ii].x, ptrSP->BC, ptrSP->sizeL);
        predNNdists[j] = dist;
    }
    prob_hit = Prob_hit(predNNdists, kill_range);
    double in_reach = 0;
    // NN which are close enough to kill selected with same prob
    in_reach = 0;
    for (j=0; j<pred->NN.size(); j++)
        if (predNNdists[j] < ptrSP->kill_range){
            prob_selected[j] = 1;
            in_reach++;
        }
    if (ptrSP->pred_kill % 5 == 2){
        int minDistIndex = std::min_element(predNNdists.begin(), predNNdists.end())
            - predNNdists.begin();
        std::fill(prob_selected.begin(), prob_selected.end(), 0);
        prob_selected[minDistIndex] = in_reach; // later divided by in_reach
    }
    if (in_reach > 0)
        vec_div221(prob_selected, in_reach);
    for (j=0; j<pred->NN.size(); j++){
        ii = pred->NN[j];
        prob_catch = prob_hit[j] * prob_selected[j] *
                     rate_catch * ptrSP->dt;
        // if(prob_catch > 0)
        //     a[ii].force[0] = 100; // marker
        a[ii].fitness -= prob_catch;
        a[ii].fit_decrease += prob_catch; // addition in case of multiple predators
        if ((ptrSP->pred_kill > 0) & (ptrSP->pred_kill < 5)){
            lucky = gsl_rng_uniform(r);
            if (prob_catch >= lucky){
                a[ii].dead = true;
                pred->kills++;
                ptrSP->Ndead++;
            }
        }
    }
}


void MoveParticle(particle &a, params * ptrSP, gsl_rng *r, double rnp){
    double h1 = 0;      // helper
    double dt = ptrSP->dt;
    double rep_strength = ptrSP->rep_strength;
    double flee_strength = ptrSP->flee_strength;
    double forcep = 0.0;
    double forcev = 0.0;
    std::vector<double> force(2, 0);
    std::vector<double> hvec(2, 0);
    double lphi = 0.0;
    int BC = ptrSP->BC;
    double sizeL = ptrSP->sizeL;
    double alg_strength = a.alg_strength;
    double speed0 = ptrSP->speed0;

    // Calc total social force
    if( a.counter_rep > 0 ){
        vec_mul221(a.force_rep, rep_strength / a.counter_rep);
        vec_mul221(a.force_alg, alg_strength / a.counter_rep);
        vec_add221(force, a.force_rep);
        vec_add221(force, a.force_alg);
    }

    if( a.counter_flee ){    // flee if pred inter
        hvec = vec_set_mag(a.force_flee, flee_strength);
        vec_add221(force, hvec);
    }
    a.force = force; // only for output
    
    // Calculate polar angle
    lphi = a.phi;

    forcep =- force[0] * sin(lphi) + force[1] * cos(lphi);

    if (a.vproj != 0)
        lphi += ( forcep * dt + rnp) / a.vproj;  // rnp = sqrt(dt * Dphi) * N(0, 1) (Wiener Process) 
    else
        lphi =  atan2(force[1], force[0]); // instantaneous direction adaption
    lphi = fmod(lphi, 2*M_PI);
    lphi = fmod(2*M_PI + lphi, 2*M_PI);   // to ensure positive angle definition
    a.phi = lphi;
    a.u[0] = cos(lphi);
    a.u[1] = sin(lphi);
    
    // Move particles with speed in units of [vel.al. range. / time]
    a.v = vec_mul(a.u, speed0);
    a.x[0] += a.v[0]*dt;
    a.x[1] += a.v[1]*dt;
    
    // Reset all forces
    // a.force[0]    =a.force[1]=0.0;
    // for output: use effective force:
    a.force_rep[0] = a.force_rep[1] = 0.0;
    a.force_alg[0] = a.force_alg[1] = 0.0;
    a.force_flee[0] = a.force_flee[1]=0.0;
    a.counter_rep = 0;
    a.counter_flee = 0;
    
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
template
void Boundary(predator &a, double sizeL,  int BC);


void MovePredator(params *ptrSP, predator *pred, std::vector<particle> &a, unsigned int s){
    // move predator based on force
    std::vector<double> com(2);
    double dist;
    std::vector<double> force(2, 0);
    double dt = ptrSP->dt;
    int N = a.size();
    unsigned int j, ii;
    std::vector<int> knn1(1);
    std::vector<double> r_pi(2);
    std::vector<int> allprey(N);
    std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
    
    // DEFINE STATE----------------------------
    std::vector<int> i_csg;     // indices of prey in circle segment
    i_csg.reserve( pred->NN.size() );
    for( j=0; j<pred->NN.size(); j++){
        ii = pred->NN[j];
        dist = CalcDist(pred->x, a[ii].x, ptrSP->BC, ptrSP->sizeL);
        if( dist <= ptrSP->predcirclerad )
            i_csg.push_back(ii);
    }
    if (pred->state == 0 && i_csg.size() > 0){
        pred->state = 1;        // change to hunt-behavior
        knn1 = GetPredKnn(a, ptrSP, pred, 1, i_csg);
        pred->cluster = GetPreyCluster(a, ptrSP, knn1[0]);
    }
    else if (pred->state == 1 && i_csg.size() == 0){
        pred->state = 0;        // change to approach-behavior
        pred->beforeCOM = true;
        if (ptrSP->pred_attack < 10)            // chooses largest cluster for approach
            pred->cluster = GetLargestCluster(a, ptrSP);
        else if (ptrSP->pred_attack < 20){      // chooses closest cluster for approach
            knn1 = GetPredKnn(a, ptrSP, pred, 1, allprey);
            pred->cluster = GetPreyCluster(a, ptrSP, knn1[0]);
        }
    }
    
    // APPROACHING STATE-----------------------
    if (pred->state == 0){
        GetCenterOfMass(a, ptrSP, pred->cluster, com);
        r_pi = CalcDistVec(pred->x, com, ptrSP->BC, ptrSP->sizeL);
        force = vec_set_mag(r_pi, 1);
    }
    
    // HUNTING STATE---------------------------
    else if (pred->state == 1){
        if ( ptrSP->pred_attack == 1 ){   // follow COM
            if (pred->beforeCOM){   // follow center of mass
                GetCenterOfMass(a, ptrSP, pred->cluster, com);
                force = CalcDistVec(pred->x, com, ptrSP->BC, ptrSP->sizeL);
                dist = vec_length(force);
                if (dist < 1 * (ptrSP->pred_speed0 + ptrSP->speed0))
                    pred->beforeCOM = false;
            }
            else  // swim straight
                force = pred->u;
        }
        else if ( ptrSP->pred_attack == 2 ){   // follow weighted com of prey with P_catch>0
            int i_prey = 0;    // index of prey with worst fitness decrease
            int min_prey = 0;  // index of closest prey 
            double mindist = N * 10 * ptrSP->rep_range;
            com[0] = com[1] = 0;  // weighted center of prey with P_catch>0
            double weight_sum = 0;
    
            for (j=0; j<i_csg.size(); j++){ // computes weighted mean position of prey in fkn (weight = 1/dist)
                ii = i_csg[j];
                if (a[ii].fit_decrease > 0){
                    com[0] += a[ii].x[0] * a[ii].fit_decrease;
                    com[1] += a[ii].x[1] * a[ii].fit_decrease;
                    weight_sum += a[ii].fit_decrease;
                }
            }
            if (weight_sum > 0)
                vec_div221(com, weight_sum);
            else{ // alternatively weight by 1/r_pi
                for (j=0; j<i_csg.size(); j++){
                    ii = i_csg[j];
                    r_pi = CalcDistVec(pred->x, a[ii].x, ptrSP->BC, ptrSP->sizeL);
                    dist = vec_length(r_pi);
                    com[0] += a[ii].x[0] / dist;
                    com[1] += a[ii].x[1] / dist;
                    weight_sum += 1 / dist;
                }
                vec_div221(com, weight_sum);
            }
            force = CalcDistVec(pred->x, com, ptrSP->BC, ptrSP->sizeL);
        }
        else           // go straight
            force = pred->u;
    }
    else{
        std::cout<< "no pred->state!" << std::endl;
        force = pred->u;
    }

    // APPLY FORCE-----------------------------------------------------
    double forcep;
    force = vec_set_mag(force, ptrSP->pred_strength);
    double lphi = pred->phi;
    double vproj = pred->vproj;     // to use correct time-step
    
    // double forcev, beta = 10;
    // forcev = force[0] * cos(lphi) + force[1] * sin(lphi);
    // pred->vproj += (beta*(speed0 - pred->vproj) + forcev) * dt;
    forcep = -force[0] * sin(lphi) + force[1] * cos(lphi);
    if (vproj != 0)
        lphi += (forcep  * dt) / vproj;
    else
        lphi =  atan2(force[1], force[1]); // instantaneous direction adaption
    lphi = fmod(lphi, 2*M_PI);
    lphi = fmod(2*M_PI + lphi, 2*M_PI);   // to ensure positive angle definition
    pred->phi = lphi;
    pred->u[0] = cos(lphi);
    pred->u[1] = sin(lphi);
    
    pred->vproj = ptrSP->pred_speed0;       // CONSTANT SPEED
    pred->v = vec_mul(pred->u, pred->vproj);
    
    // UPDATE THE POSITION---------------------------------------------
    if (!ptrSP->BC){     // only if BC=0 (periodic boundary condition)
        pred->x[0] = fmod(pred->x[0] + ptrSP->sizeL + pred->v[0] * dt, ptrSP->sizeL);
        pred->x[1] = fmod(pred->x[1] + ptrSP->sizeL + pred->v[1] * dt, ptrSP->sizeL);
    }
    else{
        pred->x[0] = pred->x[0] + pred->v[0] * dt;
        pred->x[1] = pred->x[1] + pred->v[1] * dt;
    }
    Boundary(*pred, ptrSP->sizeL, ptrSP->BC);
}


void CreatePredator(std::vector<particle> &a, params *ptrSP, predator &pred)
{
    // Creates a Predator rushing toward the particle
    std::vector<double> com(2);        // Center of mass
    unsigned int i, ii;

    // if pred already choose cluster due to memory, take this, otherwise get largest cluster
    if (pred.cluster.size() == 0){
        pred.cluster = GetLargestCluster(a, ptrSP);    // chooses largest cluster for hunt 
        ptrSP->cluster_attacked = pred.cluster;
        ptrSP->Sclu = pred.cluster.size();
    }
    GetCenterOfMass(a, ptrSP, pred.cluster, com);    // Gets the COM and puts it in com
    double avvelx = 0;                 // Average x velocity of particles
    double avvely = 0;                 // Average y velocity of particles
    for (i = 0; i < pred.cluster.size(); i++)
    {
        ii = pred.cluster[i];
        avvelx += a[ii].v[0];
        avvely += a[ii].v[1];
    }
    avvely /= pred.cluster.size();
    avvelx /= pred.cluster.size();

    // Use Average predator memory (if existend)
    if (ptrSP->pred_memo_av[0] != 0 && ptrSP->pred_memo_av[0] != 0){
        double memo_steps = (int)(ptrSP->pred_time/ptrSP->dt) + (0<fmod(ptrSP->pred_time/ptrSP->dt, 1)) -
                            (int)((ptrSP->pred_time - ptrSP->pred_memo_t)/ptrSP->dt)  - 
                            (0<fmod((ptrSP->pred_time - ptrSP->pred_memo_t)/ptrSP->dt, 1)) + 1;
        ptrSP->pred_memo_av[0] /= memo_steps;
        ptrSP->pred_memo_av[1] /= memo_steps;
        avvelx = ptrSP->pred_memo_av[0];
        avvely = ptrSP->pred_memo_av[1];
    }

    double direction = atan2(avvely, avvelx);    // Average direction of the particles in CLUSTER
    // create random angle deviation:
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    double ang_noise = ptrSP->pred_angle_noise; // normally between M_PI/6, M_PI
    std::uniform_real_distribution<> dis(-ang_noise, ang_noise);
    // TODO: side-variable only works forr ang_noise>0
    double side = dis(gen);
    if (side < 0)
        side = -1;
    else 
        side = 1;
    double angle = side * (ptrSP->pred_angle + dis(gen));  // Angle the predator is created in
    angle += direction;                // Makes angle relative to prey
    double speed = ptrSP->pred_speed0; // Predator speed
    double lphi;

    lphi = fmod(angle + M_PI, 2*M_PI);    // predator to swarm (not away)
    lphi = fmod(2*M_PI + lphi, 2*M_PI);   // to ensure positive angle definition
    pred.phi = lphi;
    pred.phi_start = lphi;     // break when pred turns back
    pred.u[0] = cos(lphi);
    pred.u[1] = sin(lphi);
    pred.v = vec_set_mag(pred.u, speed);
    pred.vproj = speed;

    //adjust predator radius to quantile distance of prey to com
    // AND save length of cluster parallel and perpend. to attack direction
    std::vector<double> r_icom(2);  // dist
    double dist;  // projected distance of r_{i,com} onto u_pred
    double maxdist = 0;  // max of projected distance of r_{i,com} onto u_pred
    for (i=0; i<pred.cluster.size(); i++){
        ii = pred.cluster[i];
        r_icom = CalcDistVec(a[ii].x, com, ptrSP->BC, ptrSP->sizeL); // a->com
        dist = r_icom[0]*pred.u[0] + r_icom[1]*pred.u[1];
        dist = fabs(dist);
        if (maxdist < dist) maxdist = dist;
    }

    // compute position of pred such that it is FOR SURE behind cluster
    // if it is in cluster the distance calculation with alpha shape does not give correct result
    pred.x[0] = com[0] - pred.u[0] *  (maxdist + 1);
    pred.x[1] = com[1] - pred.u[1] *  (maxdist + 1);
    if (!ptrSP->BC)
    {
        pred.x[0] = fmod(pred.x[0], ptrSP->sizeL);
        pred.x[1] = fmod(pred.x[1], ptrSP->sizeL);
    }
    // correct the position such that dist. of P  to alpha-shape is ptrSP->pred_radius 
    dist = DistP2AlphaShape(a, &pred, ptrSP);
    pred.x[0] += (dist - ptrSP->pred_radius) * pred.u[0];
    pred.x[1] += (dist - ptrSP->pred_radius) * pred.u[1];
    if (!ptrSP->BC)
    {
        pred.x[0] = fmod(pred.x[0], ptrSP->sizeL);
        pred.x[1] = fmod(pred.x[1], ptrSP->sizeL);
    }
}
