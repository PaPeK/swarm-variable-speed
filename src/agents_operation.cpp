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
    typedef cgPSet2::Vertex_handle                      Vertex_handle;
    typedef cgPSet2::Point                              Point;
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
    cgPSet2   PSet;      // Delauney Triangulation of Vertex set
    PSet.insert(Vr.begin(),Vr.end());
   
    cluster.push_back(id);
    unsigned int cluele = 0;        // the elements of cluster whose nn are checked
    while (cluele < cluster.size()){
      Point actual = Vr[cluster[cluele]].first;   // point of vertex whose neighbors are checked
      CGAL::Circle_2<cgK> rc(actual, squareradius);
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
    typedef cgPSet2::Vertex_handle                      Vertex_handle;
    typedef cgPSet2::Point                              Point;
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
    cgPSet2   PSet;      // Delauney Triangulation of Vertex set
    PSet.insert(Vr.begin(),Vr.end());
   
    // finds for each not sorted vertex a cluster and its other elements
    for (int i=0; i<N; i++) {
        if (vertinclu[i] < 0){
            clustertemp.resize(0);
            clustertemp.push_back(i);
            unsigned int cluele = 0;        // the elements of cluster whose nn are checked
            while (cluele < clustertemp.size()){
                Point actual = Vr[clustertemp[cluele]].first;   // point of vertex whose neighbors are checked
                CGAL::Circle_2<cgK> rc(actual, squareradius);
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


double get_elongation(std::vector<particle> &a, std::vector<double> &dir, std::vector<int> &nodes){
    std::vector<double> p_dir(2, 0); // defines perpendicular direction
    double len = vec_length(dir);
    std::vector<double> cdir = dir;
    cdir[0] /= len;
    cdir[1] /= len;
    p_dir[0] = cdir[1];
    p_dir[1] = -cdir[0];
    double min, max, p_min, p_max, dist, p_dist;
    int ii = 0;
    min = p_min = std::numeric_limits<double>::max();
    max = p_max = std::numeric_limits<double>::lowest();
    dist = p_dist = 0;
    for(int i=0; i<nodes.size(); i++){
        ii = nodes[i];
        dist = a[ii].x[0] * cdir[0] + a[ii].x[1] * cdir[1];
        p_dist = a[ii].x[0] * p_dir[0] + a[ii].x[1] * p_dir[1];
        min = fmin(min, dist);
        max = fmax(max, dist);
        p_min = fmin(p_min, p_dist);
        p_max = fmax(p_max, p_dist);
    }
    double elongation = (max - min) / (p_max - p_min);
    return fabs(elongation);
}


double AreaConvexHull(std::vector<particle> &a, std::vector<int> &nodes){
    typedef cgK::Point_2 Point_2;
    typedef CGAL::Polygon_2<cgK> Polygon_2;
    int ii;
    Point_2 points[nodes.size()];
    Point_2 convhull[nodes.size()];
    for (unsigned int i=0; i<nodes.size(); i++){
            ii = nodes[i];
            points[i] = Point_2(a[ii].x[0], a[ii].x[1]);
    }
    Point_2 *ptr = CGAL::convex_hull_2( points, points+nodes.size(), convhull);

    // create a polygon and put some points in it
    Polygon_2 p;
    for (int i=0; i<ptr-convhull; i++)
      p.push_back(convhull[i]);

    double Area = p.area();
    return Area;
}
