/*  AgentsInteract
    computes interactions between Agents based on metric or
    voronoi for SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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

#include "agents_interact.h"

void InteractionVoronoiF2F(std::vector<particle> &a, params *ptrSP)
{
    // calculates local voronoi interactions
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>         Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
    typedef Delaunay::Vertex_handle                                     Vertex_handle;
    typedef Delaunay::Edge_iterator                                     Edge_iterator;
    typedef Delaunay::Point                                             Point;
    typedef std::pair<Point, int>                                       PPoint;

    int N = a.size();

    // create the delaunay triangulation network
    std::vector< std::pair<Point, int> > Vr;   // stores the point location and index
    std::vector< std::pair< std::vector<double>, int > > posId;
    Vr.reserve(4 * N);
    posId.reserve(4 * N);
    for(int i=0; i<N; i++){
        Point p1(a[i].x[0], a[i].x[1]);
        PPoint p2 = std::make_pair(p1, i);
        Vr.push_back(p2);
        makePairAndPushBack(posId, a[i].x, i);
    }
    // produce replicate prey for periodic BC
    if (!ptrSP->BC){
        std::vector< std::pair< std::vector<double>, int > > newPosId = 
            GetCopies4PeriodicBC( posId, ptrSP->sizeL );
        for(int i=0; i<newPosId.size(); i++){
            Point p1(newPosId[i].first[0], newPosId[i].first[1]);
            PPoint p2 = std::make_pair(p1, newPosId[i].second);
            Vr.push_back(p2);
        }
    }
    // do delauney-triangulation
    Delaunay t;
    t.insert(Vr.begin(),Vr.end());

    //iterates over all finite edges and apply interaction 
    //  (infinite edges connect the infinite vertex with the vertices of the complex hull)
    for(Edge_iterator ei=t.finite_edges_begin(); ei!=t.finite_edges_end(); ei++){
      // Get a vertex from the edge, edge is stored as pair of the neighboring face and the vertex opposite to it
      Delaunay::Face& f = *(ei->first);
      int i = ei->second;
      Vertex_handle vi = f.vertex(f.cw(i));     // cw = clockwise rotation in face starting at vertex i
      Vertex_handle vj = f.vertex(f.ccw(i));    // ccw = counter clockwise .....
      IntCalcPrey(a, vi->info(), vj->info(), ptrSP, true);    // info returns the index
    }
}


void InteractionGlobal(std::vector<particle> &a, params *ptrSP)
{
    // Simple brute force algorithm for global interactions
    // checking all the N*(N-1)/2 combinations
    int i,j;
    int N = a.size();

    for(i=0;i<N;i++)
        for(j=i+1;j<N;j++)
            IntCalcPrey(a, i, j, ptrSP, true);
}


void IntCalcPrey(std::vector<particle> &a, int i, int j, params *ptrSP, bool symm)
{
    // Function updating social forces for a pair of interacting agents
    // Please Note that for local metric interaction the cutoff distance
    // is set to 1.2 x interaction range
    // check if interaction already computed (only relevant for periodic BC)
    if (!ptrSP->BC){    // BC=0: periodic BC
        int there;
        there = where_val_in_vector<unsigned int>(a[i].NN,
                                                  static_cast<unsigned int>(j));
        if (there != a[i].NN.size())   // if interaction already computed
            return;
    }
    std::vector<double> r_ji(2);
    std::vector<double> u_ji(2, 0);
    double dist_interaction;
    std::vector<double> f0(2, 0);
    std::vector<double> f1(2, 0);
    // Calc relative distance vector and corresponding unit vector
    r_ji = CalcDistVec(a[i].x, a[j].x, ptrSP->BC, ptrSP->sizeL);  // vec i->j
    dist_interaction = vec_length(r_ji);
    if(dist_interaction > 0.0)
        u_ji = vec_div(r_ji, dist_interaction);
    std::vector<double> v_ji = vec_sub(a[j].v, a[i].v);
    SFM_DRA(ptrSP, u_ji, dist_interaction, v_ji, f0, f1);
    a[i].NN.push_back(j);
    a[i].counter_rep++;
    vec_add221(a[i].force_rep, f0);
    vec_add221(a[i].force_alg, f1);
    // global, voronoi have symmetric interactions
    if (symm){
        a[j].NN.push_back(i);
        vec_mul221(f0, -1);
        vec_mul221(f1, -1);
        vec_add221(a[j].force_rep, f0);
        vec_add221(a[j].force_alg, f1);
    }
}
