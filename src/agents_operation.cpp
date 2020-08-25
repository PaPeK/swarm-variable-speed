/*  AgentsOperation
    defines operations on Agents and is part of SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#include "agents_operation.h"

void GetCenterOfMass(std::vector<particle> &a, params *ptrSP, std::vector<int> &cluster, 
                     std::vector<double> &out, bool revise, unsigned int rev_time, double quantile)
{
    // Calculates center of mass

    // if no cluster is specified, get com of all particles
    if (cluster.size() == 0){
      cluster.resize(a.size(), 0);
      std::iota (std::begin(cluster), std::end(cluster), 0); //Fill with 0, 1,...N
    }
    int NN = cluster.size();
    int i, ii;
    if (!ptrSP->BC)
    {
        // Calculates periodic boundary condition COM, see Wikipedia
        double avcosx = 0;
        double avsinx = 0;
        double avcosy = 0;
        double avsiny = 0;
        double scaling = 2*M_PI/ptrSP->sizeL;
        for (i = 0; i < NN; i++)
        {
            ii = cluster[i];
            avcosx += cos(scaling*a[ii].x[0]);
            avsiny += sin(scaling*a[ii].x[1]);
            avcosy += cos(scaling*a[ii].x[1]);
            avsinx += sin(scaling*a[ii].x[0]);
        }
        avcosx /= NN;
        avsiny /= NN;
        avcosy /= NN;
        avsinx /= NN;
        out[0] = (atan2(-avsinx,-avcosx) + M_PI)/scaling;
        out[1] = (atan2(-avsiny,-avcosy) + M_PI)/scaling;
        // Adjusts to take out outliers
        if(revise){
            std::vector<double> dist(NN);
            std::vector<double> distsort(NN);
            std::vector<double> hv(2);
            int j, i, ii = 0;
            int k;
            int counter;
            // Adjusts "rev_time" times
            for (k = 0; k < rev_time; k++)
            {
                counter = 0; // counts number of prey which are not outliers
                for (i = 0; i < NN; i++)
                    ii = cluster[i];
                    hv = CalcDistVec(a[ii].x, out, ptrSP->BC, ptrSP->sizeL);
                    dist[i] = vec_length(hv);
                distsort = dist;
                std::sort(distsort.begin(), distsort.end());
                // Makes the array dist[] contain 1 if not an outlier, 0 if it is an outlier
                for (i = 0; i < NN; i++)
                {
                    for (j = 0; j < NN; j++)
                        if (dist[i] == distsort[j])
                            break;
                    // 1 if true, 0 if false
                    dist[i] = (j < quantile*NN);
                    counter += dist[i];
                }
                // Recalculate COM
                avcosx = 0;
                avsinx = 0;
                avcosy = 0;
                avsiny = 0;
                scaling = 2*M_PI/ptrSP->sizeL;
                for (i = 0; i < NN; i++)
                {
                    if (dist[i])
                    {
                        ii = cluster[i];
                        avcosx += cos(scaling*a[ii].x[0]);
                        avsiny += sin(scaling*a[ii].x[1]);
                        avcosy += cos(scaling*a[ii].x[1]);
                        avsinx += sin(scaling*a[ii].x[0]);
                    }
                }
                avcosx /= counter;
                avsiny /= counter;
                avcosy /= counter;
                avsinx /= counter;
                out[0] = (atan2(-avsinx,-avcosx) + M_PI)/scaling;
                out[1] = (atan2(-avsiny,-avcosy) + M_PI)/scaling;
            }
        }
    }
    else
    {
        // Regular COM
        unsigned int counter = 0;
        out[1] = 0;
        out[0] = 0;
        for (i = 0; i < NN; i++)
        {
            ii = cluster[i];
            counter++;
            out[0] += a[ii].x[0];
            out[1] += a[ii].x[1];
        }
        out[0]/= counter;
        out[1]/= counter;
    }
}

std::vector<int> GetPreyCluster(std::vector<particle> &a, params *ptrSP, unsigned int id){
    // finds cluster of which particles "id" is part of, 2 particles i and j are considered connected if
    // |r_ij|<=cludist
    // Thus if they are sufficient close and sufficient parallel to each other
    typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>    Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>               Tds;
    typedef CGAL::Point_set_2<K, Tds>                 PSet2;
    typedef PSet2::Vertex_handle                      Vertex_handle;
    typedef PSet2::Point                              Point;
    typedef std::pair<Point, int>                     pp;

    double cludist = ptrSP->cludist;    // if 2 particle have dist. larger cludist -> not same cluster 
    double angdiff;                     //  angular diff. of 2 particle  
    std::vector< std::pair<Point, int> > Vr;
    int N = a.size();
    Vr.reserve(N);
    std::vector<int> vertinclu(N, -1);// for each prey the index of the 
                                      // cluster it belongs to (-1= not classified yet)
    std::vector<int> cluster;         // largest cluster
    // int icluster = 0;                 // index of largest cluster
    int clusterindex = 0;             // index of temporary cluster MUST BE POSITIVE 
    double squareradius = pow(cludist, 2);    // distance for which 2 vertexs are considered to same cluster
    unsigned int nn, there;
   
    // create triangulation point set with index
    for (int i=0; i<N; i++)
    {
          Point p1(a[i].x[0], a[i].x[1]);
          pp p2 = std::make_pair(p1, i);
          Vr.push_back(p2);
    }
    // Delaunay t;
    // t.insert(Vr.begin(),Vr.end());
    PSet2   PSet;      // Delauney Triangulation of Vertex set
    PSet.insert(Vr.begin(),Vr.end());
   
    cluster.push_back(id);
    unsigned int cluele = 0;        // the elements of cluster whose nn are checked
    while (cluele < cluster.size()){
      Point actual = Vr[cluster[cluele]].first;   // point of vertex whose neighbors are checked
      CGAL::Circle_2<K> rc(actual, squareradius);
      std::vector<Vertex_handle> NN;      // vector of nn of i
      PSet.range_search(rc, std::back_inserter(NN));
      std::vector<Vertex_handle>::const_iterator it;
      for (it=NN.begin(); it != NN.end(); it++){
        nn = (*it)->info();      // returns index of neighbor in radius
        there = where_val_in_vector<int>(cluster, nn); // returns cluster.size() if nn not found
        if (there == cluster.size() && vertinclu[nn] < 0){ 
          angdiff = fabs(a[nn].phi - a[cluster[cluele]].phi);   // ATTENTION: assuming phi in [0,2PI[
          angdiff = fabs(int(angdiff/M_PI) * M_PI - fmod(angdiff, M_PI));   // maximum angular difference is 1*Pi
          cluster.push_back(nn);
          vertinclu[nn] = clusterindex;
        }
        // UNCOMMENT FOR DEBUGGING
        // else if (there == cluster.size() && vertinclu[nn] >= 0){
        //   angdiff = fabs(a[nn].phi - a[cluster[cluele]].phi);   // ATTENTION: assuming phi in [0,2PI[
        //   angdiff = fabs(int(angdiff/M_PI) * M_PI - fmod(angdiff, M_PI));   // 
        //   if (angdiff <= cluang){
        //   std::cout << "cluster-Overlap (node(clu): " << nn << "(" <<  vertinclu[nn] << ") "
        //             << cluster[cluele] << "(" << clusterindex <<  ") phi_i-phi_j= " 
        //             <<  a[nn].phi - a[cluster[cluele]].phi << " max diff.: " << cluang << "\n";
        //   }
        // }
      }
      cluele++;
    }
    // std::cout<< "cluster.size(): " << cluster.size() << "\n";
    return cluster;
}
std::vector<int> GetLargestCluster(std::vector<particle> &a, params *ptrSP){
    // finds largest cluster from particles, 2 particles i and j are considered connected if
    // |r_ij|<=cludist
    // Thus if they are sufficient close and sufficient parallel to each other
    typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>    Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>               Tds;
    typedef CGAL::Point_set_2<K, Tds>                 PSet2;
    typedef PSet2::Vertex_handle                      Vertex_handle;
    typedef PSet2::Point                              Point;
    typedef std::pair<Point, int>                     pp;

    double cludist = ptrSP->cludist;    // if 2 particle have dist. larger cludist -> not same cluster 
    double angdiff;                     //  angular diff. of 2 particle  
    std::vector< std::pair<Point, int> > Vr;
    int N = a.size();
    Vr.reserve(N);
    std::vector<int> vertinclu(N, -1);// for each prey the index of the 
                                      // cluster it belongs to (-1= not classified yet)
    std::vector<int> cluster;         // largest cluster
    std::vector<int> clustertemp;     // cluster temporary created
    // int icluster = 0;                 // index of largest cluster
    int clusterindex = 0;             // index of temporary cluster MUST BE POSITIVE 
    double squareradius = pow(cludist, 2);    // distance for which 2 vertexs are considered to same cluster
    unsigned int nn, there;
   
    // create triangulation point set with index
    for (int i=0; i<N; i++)
    {
          Point p1(a[i].x[0], a[i].x[1]);
          pp p2 = std::make_pair(p1, i);
          Vr.push_back(p2);
    }
    // Delaunay t;
    // t.insert(Vr.begin(),Vr.end());
    PSet2   PSet;      // Delauney Triangulation of Vertex set
    PSet.insert(Vr.begin(),Vr.end());
   
    // finds for each not sorted vertex a cluster and its other elements
    for (int i=0; i<N; i++) {
        if (vertinclu[i] < 0){
            clustertemp.resize(0);
            clustertemp.push_back(i);
            unsigned int cluele = 0;        // the elements of cluster whose nn are checked
            while (cluele < clustertemp.size()){
                Point actual = Vr[clustertemp[cluele]].first;   // point of vertex whose neighbors are checked
                CGAL::Circle_2<K> rc(actual, squareradius);
                std::vector<Vertex_handle> NN;      // vector of nn of i
                PSet.range_search(rc, std::back_inserter(NN));
                std::vector<Vertex_handle>::const_iterator it;
                for (it=NN.begin(); it != NN.end(); it++){
                    nn = (*it)->info();      // returns index of neighbor in radius
                    there = where_val_in_vector<int>(clustertemp, nn); // returns clustertemp.size() if nn not found
                    if (there == clustertemp.size() && vertinclu[nn] < 0){ 
                        angdiff = fabs(a[nn].phi - a[clustertemp[cluele]].phi);   // ATTENTION: assuming phi in [0,2PI[
                        angdiff = fabs(int(angdiff/M_PI) * M_PI - fmod(angdiff, M_PI));   // maximum angular difference is 1*Pi
                        clustertemp.push_back(nn);
                        vertinclu[nn] = clusterindex;
                    }
                    // UNCOMMENT FOR DEBUGGING
                    // else if (there == clustertemp.size() && vertinclu[nn] >= 0){
                    //   angdiff = fabs(a[nn].phi - a[clustertemp[cluele]].phi);   // ATTENTION: assuming phi in [0,2PI[
                    //   angdiff = fabs(int(angdiff/M_PI) * M_PI - fmod(angdiff, M_PI));   // 
                    //   if (angdiff <= cluang){
                    //   std::cout << "cluster-Overlap (node(clu): " << nn << "(" <<  vertinclu[nn] << ") "
                    //             << clustertemp[cluele] << "(" << clusterindex <<  ") phi_i-phi_j= " 
                    //             <<  a[nn].phi - a[clustertemp[cluele]].phi << " max diff.: " << cluang << "\n";
                    //   }
                    // }
                }
                cluele++;
            }
            if (clustertemp.size() > cluster.size()){
              cluster.resize(0);
              cluster = clustertemp;
              // icluster = clusterindex;
            }
            clusterindex++;   // give next cluster another index
            }
    }
    // ptrSP->Nclu = clusterindex; 
    // ptrSP->NLclu = cluster.size();
    // std::cout<< "cluster.size(): " << cluster.size() << "\n";
    return cluster;
}

std::vector<unsigned int> GetPredVoronoiNN(std::vector<particle> &a, 
                                           params *ptrSP, predator *pred,
                                           std::vector<int> &nodes)
{
    // voronoi nearest neighbors of predator, only considers prey whose index is in "nodes"
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>         Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
    typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
    typedef Delaunay::Vertex_handle                                     Vertex_handle;
    typedef Delaunay::Point                                             Point;
    typedef std::pair<Point, int>                                       PPoint;
    std::vector<unsigned int> voronoinn; // contains indices of prey in circle(r= circlerad) around pred
    voronoinn.reserve(7);
    std::vector< std::pair<Point, int> > Vr;   // stores the point location and index

    int N = nodes.size();       // # of prey considered
    Vr.reserve(N+1);
    for (int i=0; i<N; i++){
        Point p1(a[nodes[i]].x[0], a[nodes[i]].x[1]);
        PPoint p2 = std::make_pair(p1, nodes[i]);  // prey labeled with corresponding index
        Vr.push_back(p2);
    }

    Point p1(pred->x[0], pred->x[1]);
    PPoint p2 = std::make_pair(p1, -1);     // predator labeled with negative index
    Vr.push_back(p2);

    Delaunay t;
    t.insert(Vr.begin(),Vr.end());

    Point actual = Vr[Vr.size()-1].first;
    Vertex_handle v =  t.nearest_vertex(actual);
    Vertex_circulator vc = t.incident_vertices(v),
    done(vc);
    
    if (vc != 0) {
      do{ 
        if (!t.is_infinite(vc)){    // avoids the infinite vertex (marks convex hull nodes)
        voronoinn.push_back(vc->info());
        }
      }while(++vc != done);
    }
    return voronoinn;
}


std::vector<int> GetPredKnn(std::vector<particle> &a, params *ptrSP, predator *pred, unsigned int K, std::vector<int> &nodes)
{
    // returns vector of indices to K nearest neigbors of predator, 
    // if inclu: only to prey in cluster
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_2                                     Point_2;
    typedef boost::tuple<Point_2,int>                           Point_and_int;
    typedef CGAL::Search_traits_2<Kernel>                       Traits_base;
    typedef CGAL::Search_traits_adapter<Point_and_int,
      CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
      Traits_base>                                              Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
    typedef K_neighbor_search::Tree                             Tree;

    // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
    std::vector<Point_2> points;
    std::vector<int>     indices;
    int partner;
    K += 1;    // because the cgal also lists the node itself
    int N = nodes.size();       // # of prey considered
    if ( K > N )
        K = N + 1;  // +1 because of the predator

    for (int i=0; i<N; i++){
        Point_2 p1(a[nodes[i]].x[0], a[nodes[i]].x[1]);
        points.push_back(p1);
        indices.push_back(nodes[i]);       // prey has positive index
    }
    
    std::vector<int>  results(K-1);   // save index of the K nearest neighbors
    Point_2 p2(pred->x[0], pred->x[1]);
    points.push_back(p2);
    indices.push_back(-1);              // predator has negative index

    Tree tree(
      boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
      boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )  
    );

    Point_2 query(pred->x[0], pred->x[1]);  // defines point in 2D-space from which NN shall be found
    // search K nearest neighbours
    K_neighbor_search search(tree, query, K);
    int i = 0;
    for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
    {
      // it->first: tuple of (point, index), it->second:(distance to query)^2
      partner = boost::get<1>(it->first);     // acces first element which is tuple of (point, index)-> get<1>=index
      if (partner >= 0){
              results[i] = partner;
              i++;
      }
    }
    return results;
}


std::vector< CGAL::Exact_predicates_inexact_constructions_kernel::Segment_2 >
        AlphaShapeSegments(std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Point_2> &points,
                           double r){
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::FT FT;
    typedef K::Point_2  Point;
    typedef K::Segment_2  Segment;
    typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
    typedef CGAL::Alpha_shape_face_base_2<K>  Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;
    typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
    typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
    Alpha_shape_2 A(points.begin(), points.end(),
                    FT(r * r),      // alpha-value=sqrt(r)
                    Alpha_shape_2::GENERAL); // REGULARIZED only include edges which are part of a 2D-face
    std::vector<Segment> segments;
    // alpha_edges( A, std::back_inserter(segments));  // UNORDERED segments of alpha-shape
    for(Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin();
        it!=A.alpha_shape_edges_end();
        ++it){
        segments.push_back(A.segment(*it));
    }
    return segments;
}


std::vector<double> Dist2AlphaShape(std::vector<particle> &a,
                        params *ptrSP){
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2  Point;
    typedef K::Segment_2  Segment;
    std::vector<Point> points;
    points.reserve(a.size());
    double dist;
    // compute the min and max values for ray-space
    for (unsigned int i=0; i<a.size(); i++)
        points.push_back(Point(a[i].x[0], a[i].x[1]));
    double r = ptrSP->cludist/2. + 0.001; // radius of ice-cone (must hold 2*r > cludist!!!)
    std::vector<Segment> segments = AlphaShapeSegments(points, r);
    // rc::vec_t x0(segments[i][0][0], segments[i][0][1]);
    std::vector<double> x(2);
    std::vector<double> distances(a.size());
    double min_dist = 0;
    for (unsigned int i=0; i<a.size(); i++){
        x = a[i].x;
        min_dist = a.size() * a.size(); // initialize hughe distance
        for (unsigned int j=0; j<segments.size(); j++){
            std::vector<double> A {segments[j][0][0], segments[j][0][1]};
            std::vector<double> B {segments[j][1][0], segments[j][1][1]};
            dist = DistPoint2Segment(x, A, B);
            if (dist < min_dist)
                min_dist = dist;
            if (dist == 0)
                break;
        }
        distances[i] = min_dist;
    }
    return distances;
}

double DistP2AlphaShape(std::vector<particle> &a, predator *pred,
                        params *ptrSP){
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2  Point;
    typedef K::Segment_2  Segment;
    std::vector<Point> points;
    points.reserve(pred->cluster.size());
    int ii;
    std::vector<rc::line_t> alphaShape_lines;
    double dist;
    double epsilon = 0.1; // ensures that no point lies on the boundary IMPORTANT!!!!!!!!!!!!!!!!!!
    // compute the min and max values for ray-space
    double min_x, max_x, min_y, max_y;
    min_x = max_x = pred->x[0];
    min_y = max_y = pred->x[1];
    for (unsigned int i=0; i<pred->cluster.size(); i++){
        ii = pred->cluster[i];
        if (min_x > a[ii].x[0]) min_x = a[ii].x[0];
        if (max_x < a[ii].x[0]) max_x = a[ii].x[0];
        if (min_y > a[ii].x[1]) min_y = a[ii].x[1];
        if (max_y < a[ii].x[1]) max_y = a[ii].x[1];
    }
    double offset[2] = {min_x - epsilon, min_y - epsilon};
    for (unsigned int i=0; i<pred->cluster.size(); i++){
            ii = pred->cluster[i];
            points.push_back(Point(a[ii].x[0] - offset[0], a[ii].x[1] - offset[1]));
    }
    double r = ptrSP->cludist/2. + 0.001; // radius of ice-cone (must hold 2*r > cludist!!!)
    std::vector<Segment> segments = AlphaShapeSegments(points, r);
    
    // std::ofstream outfile("temp.dat");
    // create lines from convex hull:
    for (unsigned int i=0; i<segments.size(); i++){
        // outfile<<  segments[i] << "\n";
        rc::vec_t x0(segments[i][0][0], segments[i][0][1]);
        rc::vec_t x1(segments[i][1][0], segments[i][1][1]);
        alphaShape_lines.push_back(rc::line_t(x0, x1));
    }

    // create size of ray-space:
    // (always 16 bins in each dimension, regardless of the size)
    rc::vec_t size(max_x - min_x + 2 * epsilon,
                   max_y - min_y + 2 * epsilon);
    // create ray-space:
    fov::ray_space rspace(size);
    // create bodies and add them to ray-space
    rc::body_t *body = rspace.create_body_from_sequence(alphaShape_lines, 0);
    rspace.add_body_to_space(body);

    // create container for bodies found by raycasting
    rc::body_t ** visible_bodies;
	visible_bodies  = new rc::body_t*[1];

    // define ray start and end points:
	double kMaxLength = sqrt(size.x*size.x + size.y*size.y);
    rc::vec_t orig(pred->x[0] - offset[0],
                   pred->x[1] - offset[1]);
    rc::vec_t end(pred->x[0] - offset[0] + kMaxLength * pred->u[0],
                  pred->x[1] - offset[1] + kMaxLength * pred->u[1]);

    rc::vec_t terminus = rspace.terminus_of_ray(orig, end, 0.0, &visible_bodies[0]);
    // outfile<< orig.x << " " << orig.y << " " <<  terminus.x << " " << terminus.y << "\n";
    // outfile<< orig.x << " " << orig.y << " " <<  end.x << " " << end.y << "\n";
    rc::vec_t distvec(terminus, orig);
    dist = distvec.magnitude();

    rspace.reset();
    delete [] visible_bodies;

    return dist;
}


template <class O, class I>
std::vector<O> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<I> &nodes)
{
    // returns vector of indices of prey in fron of pred 
    // only prey considered whose index is in "nodes" 
    unsigned int i;
    I ii;
    std::vector<double> r_pi(2);
    double front;   // if positive -> prey is in front 
    std::vector<O> results;
    for(i=0; i<nodes.size(); i++){
        ii = nodes[i];
        r_pi = CalcDistVec(pred->x, a[ii].x, ptrSP->BC, ptrSP->sizeL);
        front = pred->u[0] * r_pi[0] + pred->u[1] * r_pi[1];
        if (front > 0.0)
            results.push_back(ii);
    }
    return results;
}
template
std::vector<int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<int> &nodes);
template
std::vector<unsigned int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<unsigned int> &nodes);
template
std::vector<int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<unsigned int> &nodes);
template
std::vector<unsigned int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<int> &nodes);


void split_notInCluster(std::vector<particle> &a, std::vector<particle> &d,
                        std::vector<int> &cluster, std::vector<predator> &preds){
    int removed = 0;
    int ii = 0; 
    std::sort(cluster.begin(), cluster.end()); // sorted in ascending order
    std::reverse(cluster.begin(), cluster.end()); // sorted in descending order
    int i_cluster = 0;
    std::vector<int> pred_cluster = preds[0].cluster;
    for (int i=a.size()-1; i>=0; i--){
        // not in cluster -> remove
        // OR
        // whole cluster already checked -> not in cluster -> remove
        if (i != cluster[i_cluster] || i_cluster >= cluster.size()){
            // for debugging
            if ((i < cluster[i_cluster]) && (i_cluster < cluster.size()))
                std::cout<< "ATTENTION!!! THIS IS WRONG!!" << 
                    a.size() << ", " << d.size() << ", " << cluster.size()<< 
                    ",, " << i << ", " << i_cluster << ", " << cluster[i_cluster] << std::endl;
            d.push_back(a[i]);
            a.erase(a.begin() + i);
            auto found2 = std::find(pred_cluster.begin(), pred_cluster.end(), i);
            if(pred_cluster.end() != found2)
                pred_cluster.erase(found2);
            // correct indices of existing cluster 
            for (int j=0; j<pred_cluster.size(); j++){
                if (pred_cluster[j] > i)
                    pred_cluster[j]--;
            }
        }
        else    // if it is in cluster
            i_cluster++;
    }
    for (unsigned int i=0; i<preds.size(); i++)
        preds[i].cluster = pred_cluster;
}


void split_dead(std::vector<particle> &a, std::vector<particle> &d,
                std::vector<predator> &preds){
    for (int i=0; i<a.size(); i++){
        if (a[i].dead){
            d.push_back(a[i]);
            a.erase(a.begin() + i);
            for (int j=0; j<preds.size(); j++){
                auto found = std::find(preds[j].cluster.begin(), preds[j].cluster.end(), i);
                if(preds[j].cluster.end() != found)
                    preds[j].cluster.erase(found);
                // correct indices of existing cluster 
                for (int k=0; k<preds[j].cluster.size(); k++){
                    if (preds[j].cluster[k] > i)
                        preds[j].cluster[k] -= 1;
                }
            }
            i--;
        }
    }
}

void merge_dead(std::vector<particle> &a,
                                 std::vector<particle> &d){
    std::vector<particle> all(a.size() + d.size());
    for (int i=0; i<a.size(); i++)
        all[a[i].id] = a[i];
    for (int i=0; i<d.size(); i++)
        all[d[i].id] = d[i];
    d.resize(0);
    a = all;
}


void makePairAndPushBack(std::vector< std::pair< std::vector<double>, int > > &vecpair,
                         std::vector<double> &vec, int id){
    std::pair< std::vector<double>, int > newPair;
    newPair = std::make_pair( vec, id );
    vecpair.push_back( newPair );
}


std::vector< std::pair< std::vector<double>, int > > GetCopies4PeriodicBC(
        std::vector< std::pair< std::vector<double>, int > > &posId, double L)
{
    bool lower;
    bool left;
    double x, y;
    int id;
    std::vector<double> pos;
    std::vector<double> right {L, 0};
    std::vector<double> up {0, L};
    std::vector<double> newPos;
    std::vector< std::pair< std::vector<double>, int > > newPosId;
    newPosId.reserve(3 * posId.size()); // shifted copies of original pos (SAME ID)
    for (int i=0; i < posId.size(); i++){
        pos = posId[i].first;
        id = posId[i].second;
        if (pos[0] < L/2)
            left = true;
        if (pos[1] < L/2)
            lower = true;
        if (left && lower){
            newPos = vec_add(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_add221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_add(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else if (left && !lower){
            newPos = vec_sub(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_add221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_add(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else if (!left && !lower){
            newPos = vec_sub(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_sub221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_sub(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else{ // (!left && lower)
            newPos = vec_add(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_sub221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_sub(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
    }
    return newPosId;
}
