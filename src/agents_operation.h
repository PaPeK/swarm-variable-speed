/*  AgentsOperation
    defines operations on Agents. 
    This code is part of SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#ifndef agents_operation_H
#define agents_operation_H
// OWN MODULES:
#include "agents.h"
#include "mathtools.h"
#include "settings.h"
// for raycasting (slightly modified code from Colin Twomey)
#include "ray_casting.h"
#include "geometry.h"
#include "fov.h"

#include <set>
#include <limits>
#include <algorithm>    // to use std::set_difference(...)
#include <cassert>      // to use assert
// for nearest neighbor search needed
#include <utility>
// for CGAL functionality:
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>      // Trangulation Class with member incident_vertices
#include <CGAL/convex_hull_2.h>     // computes convex hull
#include <CGAL/Polygon_2.h>         // compute area of convex hull, of voro-cells
#include <CGAL/Point_set_2.h>                   // Trangulation Class with member nearest_neighbors
#include <CGAL/Triangulation_2.h>               // for Trangulation needed
#include <CGAL/Triangulation_vertex_base_with_info_2.h> // for Trangulation needed
#include <CGAL/Search_traits_2.h>               // for nearest neighbor search needed
#include <CGAL/Search_traits_adapter.h>         // for nearest neighbor search to change points class 
#include <CGAL/Orthogonal_k_neighbor_search.h>  // for nearest neighbor search needed
#include <CGAL/property_map.h>                  // for nearest neighbor search needed
#include <boost/iterator/zip_iterator.hpp>      // for nearest neighbor search needed

void GetCenterOfMass(std::vector<particle> &, params *,
                     std::vector<int> &, std::vector<double> &out,
                     bool revise=false, unsigned int rev_time=1, double quantile=0.9); // gives the center of mass of cluster
std::vector<int> GetLargestCluster(std::vector<particle> &, params *); //finds indices of largest cluster
std::vector<int> GetPreyCluster(std::vector<particle> &a, params *ptrSP, unsigned int id); //finds indices of cluster where prey "id" is in
std::vector< CGAL::Exact_predicates_inexact_constructions_kernel::Segment_2 >
        AlphaShapeSegments(std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Point_2> &points,
                           double r);
std::vector<double> Dist2AlphaShape(std::vector<particle> &a,
                        params *ptrSP);
// returns indicese of Prey(in "nodes") in front of Pred
void makePairAndPushBack(std::vector< std::pair< std::vector<double>, int > > &vecpair,
                         std::vector<double> &vec, int id);
std::vector< std::pair< std::vector<double>, int > > GetCopies4PeriodicBC(
        std::vector< std::pair< std::vector<double>, int > > &posId, double L);
#endif
