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
#ifndef agents_interact_H
#define agents_interact_H
// OWN MODULES:
#include "agents.h"
#include "agents_operation.h"
#include "mathtools.h"
#include "social_forces.h"
#include "settings.h"

#include <algorithm>
#include <vector>

// for nearest neighbor search needed
#include <utility>
// for raycasting (slightly modified code from Colin Twomey)
#include "ray_casting.h"
#include "geometry.h"
#include "fov.h"
// for CGAL functionality:
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>      // Trangulation Class with member incident_vertices
#include <CGAL/convex_hull_2.h>     // computes convex hull
#include <CGAL/Polygon_2.h>         // compute area of convex hull
#include <CGAL/Point_set_2.h>                   // Trangulation Class with member nearest_neighbors
#include <CGAL/Triangulation_2.h>               // for Trangulation needed
#include <CGAL/Triangulation_vertex_base_with_info_2.h> // for Trangulation needed
#include <CGAL/Search_traits_2.h>               // for nearest neighbor search needed
#include <CGAL/Search_traits_adapter.h>         // for nearest neighbor search to change points class 
#include <CGAL/Orthogonal_k_neighbor_search.h>  // for nearest neighbor search needed
#include <CGAL/property_map.h>                  // for nearest neighbor search needed
#include <boost/iterator/zip_iterator.hpp>      // for nearest neighbor search needed

// INTERACTION COMPUTATION-------------------------------------------
// voronoi: fish-fish
void InteractionVoronoiF2F(std::vector<particle> &a, params *);
// voronoi: fish-fish, fish-pred
void InteractionVoronoiF2FP(std::vector<particle> &, params *,
        std::vector<predator> &, bool dummy = false);
// global: fish-fish
void InteractionGlobal(std::vector<particle> &, params *);
// global: fish-fish, fish-pred
void InteractionPredGlobal(std::vector<particle> &, params *,
        std::vector<predator> &);
// fish-fish:
void IntCalcPrey(std::vector<particle> &, int, int, params *, bool symm);
// fish-pred:
void IntCalcPred(std::vector<particle> &, int, predator &, params *);

#endif
