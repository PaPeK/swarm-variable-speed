/*
 *  ray_casting.h
 *  This file is part of FOVEA
 *
 *  Created by Colin Twomey on 02.03.2012
 *  Copyright (c) 2012 Colin Twomey
 *
 *  FOVEA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FOVEA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with FOVEA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FOV_RAY_CASTING_H
#define FOV_RAY_CASTING_H

#include <set>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream> // HELPER

#include "geometry.h"
#include "spatial_hash.h"

namespace rc {
using namespace gm;

typedef struct {
	bool   hit;
	vec_t  hit_point;
	vec_t  hit_normal;  // PPK: normal vector to surface where hit appeared
	line_t *line;       // points to line (contains info regarding body)
} ray_result_t;


/*
 *  calc_outcode and collide_linesegment_boundingbox code adapted from
 *  Schroeder, Tim, "Collision Detection Using Ray Casting,"
 *  Game Developer Magazine, pp. 50-57, August 2001
 */

// Cohen-Sutherland clipping outcodes
#define CLIP_RIGHT		(1<<0)
#define CLIP_LEFT		(1<<1)
#define CLIP_TOP		(1<<2)
#define CLIP_BOTTOM		(1<<3)

bool collide_linesegment_boundingbox(vec_t &bbox_min, vec_t &bbox_max,
                                     vec_t &p1, vec_t &p2, vec_t &intercept);
// static ulong calc_outcode(vec_t &bbox_min, vec_t &bbox_max, vec_t &pnt, // PPK: commented out "static" because it was not used (compiler-warning)
ulong calc_outcode(vec_t &bbox_min, vec_t &bbox_max, vec_t &pnt,
                          double margin=0.0);


/*
 *  segment_intersect, line_intersect, and classify_point code adapted
 *  from Laszlo, Michael J., "Computational Geometry and Computer
 *  Graphics in C++" (1996)
 */
typedef enum {
	LEFT, RIGHT, BEYOND, BEHIND, BETWEEN, ORIGIN, DESTINATION
} point_pos_type_t;

typedef enum {
	COLLINEAR, PARALLEL, SKEW, SKEW_CROSS, SKEW_NO_CROSS, SKEW_NO_CROSS_2
} cross_type_t;

cross_type_t     segment_intersect (line_t &v, line_t &u, vec_t *x);
cross_type_t     line_intersect    (line_t &v, line_t &u, double &s);
point_pos_type_t classify_point    (vec_t &pt, line_t &edge);


/*
 *  A 2d spatial hashing class, with an efficient quadtree
 *  overlay for ray-casting.
 */
class space : public spatial_hash<2, line_t*> {
public:
	space(SHuint dim[2], SHfloat size[2], SHfloat offset[2])
		: spatial_hash<2, line_t*>(dim, size, offset),
		  m_body_counter(0)
	{ }
	
	body_t* add_body     (uint d, uint owner_tid=-1);   // PPK: creates body_t pointer heaped body_t and creates the heaped-arrays of line_t which define the body  and stores the body in m_bodes (std::list<body_t*>)
	void    attach_body  (body_t *body);
	void    clear_bodies (void);    // PPK: goes through all body_t in m_body and deletes the heaped-array of line_t and the heaped body_t

	ray_result_t cast_ray(line_t ray);
	ray_result_t cast_ray_from_body(line_t ray, uint body_id);
	ray_result_t cast_ray_neglect_bodies(line_t ray,
                                         std::vector<uint> &bodies_id);
	void    bins_between(SHfloat a[2], SHfloat b[2], SHuint h_b,
	                     std::set<SHuint> *result);

private:
	std::list<body_t*> m_bodies;
	uint m_body_counter;

	double traverse_ray(vec_t &a, vec_t end, uint h_a);
};


}; // rc

#endif // FOV_RAY_SPACE_H
