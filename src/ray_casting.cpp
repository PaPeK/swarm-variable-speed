/*
 *  ray_casting.cpp
 *  This file is part of FOVEA
 *
 *  Created by Colin Twomey on 02.03.2011
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

#include <float.h>
#include <set>
#include <utility>
#include <iostream>

#include "ray_casting.h"


namespace rc {

const double kErrorMargin     = 0.1;
const double kTraversalOffset = 0.0001;

/*
 *  calculates the cohen-sutherland outcode for a point and a bounding box.
 *  
 *  bbox_min:	min vector of the bounding box
 *  bbox_max:	max vector of the bounding box
 *  pnt:		the point to check
 *  
 *  returns:	the outcode
 *          CLIP_RIGHT, CLIP_LEFT as expected
 *          CLIP_TOP, CLIP_BOTTOM should be exchanged....???
 */
// static unsigned long calc_outcode(   // PPK: commented out "static" because it was not used (compiler-warning)
unsigned long calc_outcode(
	vec_t &bbox_min,
	vec_t &bbox_max,
	vec_t &pnt,
	double margin)
{
	ulong outcode = 0;  // bitwise 0000

	if (pnt.x > bbox_max.x+margin)
		outcode |= CLIP_RIGHT;  // | is bitwise OR: a = 100 |= 110 -> a = a | 110 = 110 
	else if (pnt.x < bbox_min.x-margin)
		outcode |= CLIP_LEFT;

	if (pnt.y > bbox_max.y+margin)
		outcode |= CLIP_BOTTOM;
	else if (pnt.y < bbox_min.y-margin)
		outcode |= CLIP_TOP;

	return outcode;
}

/*
 *  Determines if a linesegment intersects a bounding box. this is based on
 *  the Cohen-Sutherland line-clipping algorithm.
 *  
 *   bbox_min:  bounding box min vector
 *   bbox_max:  bounding box max vector
 *         p1:  end point of line segment
 *         p2:  other end point
 *  intercept:  (out) the point in/on the bounding box where the intersection 
 *              occured.  Note that this point may not be on the surface
 *              of the box.
 *  
 *    returns:  true if the segment and box intersect.
 */
bool collide_linesegment_boundingbox(
	vec_t &bbox_min,
	vec_t &bbox_max,
	vec_t &p1,
	vec_t &p2,
	vec_t &intercept)
{
	ulong outcode1,
	      outcode2;

	outcode1 = calc_outcode(bbox_min, bbox_max, p1);
	outcode2 = calc_outcode(bbox_min, bbox_max, p2);

	if (outcode1 == 0 && outcode2 == 0) {
		// line segment inside bounding box
		intercept = p1;
		return false;
	}

	if ((outcode1 & outcode2) > 0 ) {
		// both points on same side of box
		return false;
	}

	// check intersections
    // PPK: here only checked if p1 is outside the box
    //         the same code below needs to be repeated for p2 
    //         (probably not done because code not in use)
	if (outcode1 & (CLIP_RIGHT | CLIP_LEFT)) {
		if (outcode1 & CLIP_RIGHT)
			intercept.x = bbox_max.x;
		else
			intercept.x = bbox_min.x;

		float x1 = p2.x - p1.x;
		float x2 = intercept.x - p1.x;
		intercept.y = p1.y + x2 * (p2.y - p1.y) / x1;

		if (intercept.y <= bbox_max.y && intercept.y >= bbox_min.y)
			return true;
	}

	if (outcode1 & (CLIP_TOP | CLIP_BOTTOM)) {
		if (outcode1 & CLIP_TOP)
			intercept.y = bbox_max.y;
		else
			intercept.y = bbox_min.y;

		float y1 = p2.y - p1.y;
		float y2 = intercept.y - p1.y;
		intercept.x = p1.x + y2 * (p2.x - p1.x) / y1;

		if (intercept.x <= bbox_max.x && intercept.x >= bbox_min.x)
			return true;
	}

	// nothing found
	return false;
}

/*PPK:
 *classifies the position of a point relative to a line
 * OUTPUT:  LEFT, RIGHT
 *          BEHIND  point before start of line
 *          BEYOND  point behind end of line
 *          ORIGIN  point at start of line
 *          DESTINATION  point at end of line
 */
point_pos_type_t classify_point(rc::vec_t &pt, rc::line_t &edge)
{
// PPK: difference vector of V(=edge): a = (\Delta Edge_x, \Delta Edge_y)
	vec_t a = vec_t(edge.end.x - edge.start.x,
	                edge.end.y - edge.start.y);
    // PPK: difference vector of vector U with: start = edge_start, end = pt
	vec_t b = vec_t(pt.x - edge.start.x,
	                pt.y - edge.start.y);

    // PPK: comparing the slopes: sa is numerator of (m_U - m_V)
	double sa = a.x * b.y - b.x * a.y;

	if (sa > 0.0) return LEFT;  // PPK: m_U slope is larger -> left
	if (sa < 0.0) return RIGHT;  // PPK: m_U slope is smaller -> right 
	if ((a.x * b.x < 0.0) || (a.y * b.y < 0.0))
		return BEHIND;
	if (a.magnitude() < b.magnitude())
		return BEYOND;
	if (edge.start.x == pt.x && edge.start.y == pt.y)
		return ORIGIN;
	if (edge.end.x == pt.x && edge.end.y == pt.y)
		return DESTINATION;
	return BETWEEN;
}


/*PPK:
 *given two lines v, u it computes if they 
 *  PARALLEL, COLLINEAR
 *  SKEW and computes 
 *  s which determines if the two lines do cross or not
 */
cross_type_t line_intersect(line_t &v, line_t &u, double &s)
{
	vec_t a = v.start; vec_t b = v.end;
	vec_t c = u.start; vec_t d = u.end;
    // PPK: n is 'strange' difference vector of u: n = (\Delta u_y, -\Delta u_x)
	vec_t n = vec_t(d.y - c.y, c.x - d.x);

    // PPK: denom is numerator of m_u-m_v
	double denom = n.x * (b.x - a.x) + n.y * (b.y - a.y);

	if (denom == 0.0) {
		point_pos_type_t ppt = classify_point(v.start, u);
		if (ppt == LEFT || ppt == RIGHT)
			return PARALLEL;
		return COLLINEAR;
	}

    // PPK: num is numerator of m_u-m_w with vector w connecting both start points
	double num = n.x * (a.x - c.x) + n.y * (a.y - c.y);
    // PPK: no idea of mathmatical derivation but if s<0 or s>1 -> no cross
    //                                               0 <= s <=1 -> cross
	s = -num / denom;

	return SKEW;
}


cross_type_t segment_intersect(line_t &v, line_t &u, vec_t *x)
{
	double t;
	cross_type_t ct = line_intersect(u, v, t);
	if (ct == COLLINEAR || ct == PARALLEL)
		return ct;
	if (t < 0.0 || t > 1.0)
		return SKEW_NO_CROSS;

    // PPK: don't know why recompute ct and exchanging u<->v
    //      ... did not understand what t is
	ct = line_intersect(v, u, t);
	if (0.0 <= t && t <= 1.0) {
     	vec_t vn = vec_t(v.end.x - v.start.x, v.end.y - v.start.y);
     	x->x = v.start.x + vn.x * t;
     	x->y = v.start.y + vn.y * t;

		return SKEW_CROSS;
	}

	return SKEW_NO_CROSS_2;
}


body_t* space::add_body(uint d, uint owner_tid)
{
	body_t *body = new body_t;
	m_bodies.push_back(body);

	body->owner_tid   = owner_tid;
	body->total_lines = d;
	body->lines       = new line_t[d];
	body->bid         = m_body_counter;
	m_body_counter++;

	return body;
}


void space::clear_bodies(void)
{
     while (!m_bodies.empty()) {
        body_t *body = m_bodies.back();
        m_bodies.pop_back();

		if (body != NULL) {
			if (body->lines != NULL)
	        	delete [] body->lines;
	        delete body;
		}
     }
}


void space::attach_body(body_t *body)
{
	for (uint i = 0; i < body->total_lines; i++) {
     	SHfloat a[] = { body->lines[i].start.x, body->lines[i].start.y };
     	SHfloat b[] = { body->lines[i].end.x,   body->lines[i].end.y   };

     	std::set<SHuint> hs;
        // PPK: saves all bins between a and b in hs
     	bins_between(a, b, hash(b), &hs);

        // PPK: places the line in all the bins it crosses
		for (std::set<SHuint>::iterator j = hs.begin();
			 j != hs.end(); j++)
		{
			place(*j, body->lines + i);
		}
	}
}


/*PPK:
 * cast ray by checking 
 *  if any line in box of start point intersects with ray:
 *      if yes: check if intersection is in box
 *          if yes: mark hit BUT check all other lines as well
 *                  to find shortest hit
 *  if no: move start point to next box till end of ray point
 */
ray_result_t space::cast_ray(line_t ray)
{
	vec_t &start = ray.start;
	vec_t &end   = ray.end;
	
	ray_result_t result;
	result.hit_point = vec_t::kZero;
	result.hit       = false;
	result.line      = NULL;

	SHfloat s[] = { start.x, start.y };
	SHfloat e[] = { end.x,     end.y };

	SHfloat moved_by = 0.0;
	SHuint h_s;

	// keep track of bins traversed
	std::set<SHuint> traversed;

	do {
		h_s = hash(s);

		// check for any strange traversal patterns
		if (h_s == SH_OUT_OF_BOUNDS) break;
		else if (traversed.count(h_s) > 0) break;
		else traversed.insert(h_s);

		std::set<line_t*> *lines = index(h_s);
		
		if (!lines->empty()) {
			double nearest = DBL_MAX;
			for (std::set<line_t*>::iterator i = lines->begin();
			     i != lines->end(); i++)
			{
				vec_t hit_point = vec_t::kZero;
				line_t   r_copy = ray;
				line_t   c_copy = **i;
				cross_type_t ct =
					segment_intersect(r_copy, c_copy, &hit_point);

				// check for an intersection
				if (ct == SKEW_CROSS) {
					SHfloat min_pt[2] = {};
					SHfloat max_pt[2] = {};
					index_bounds(h_s, min_pt, max_pt);
					vec_t min_v(min_pt[0], min_pt[1]),
					      max_v(max_pt[0], max_pt[1]);

					// check that the intersection is in the box
					if (calc_outcode(min_v, max_v, hit_point, kErrorMargin)
					    == 0)
					{
						double dist = vec_t::distance(ray.start, hit_point);
						if (dist < nearest) {
							nearest          = dist;
							result.hit_point = hit_point;
							result.hit       = true;
							result.line      = *i;

							point_pos_type_t pt 
								= classify_point(start, c_copy);

                            // PPK: compute normal vector to colliding line 
							if (pt == LEFT) {
								result.hit_normal = vec_t(
									c_copy.start.y - c_copy.end.y,
									c_copy.end.x   - c_copy.start.x
								);
							} else if (pt == RIGHT) {
								result.hit_normal = vec_t(
									c_copy.end.y   - c_copy.start.y,
									c_copy.start.x - c_copy.end.x
								);
							}

							result.hit_normal.normalize();
						}
					}
				}
			}
		}

		if (result.hit != true) {
			vec_t current(s[0], s[1]);
			vec_t end(e[0], e[1]);

            // PPK: moved_by==0 ONLY if "current" and "end" are same point 
			moved_by = traverse_ray(current, end, h_s);
			s[0] = current.x;
			s[1] = current.y;
		}
	} while (moved_by > 0.0 && result.hit != true);

	return result;
}


/*PPK:
 * computes the interception of line a->end and the hash-box
 * which contains a(= the start point)
 *
 * OUTPUT:
 *  double  distance between a and interception
 *  a       modifies start point a such that it is now in new bin
 *              after the interception with the spatial-box
 */
double space::traverse_ray(vec_t &a, vec_t end, uint h_a)
{
 	SHfloat min_pt[2] = {};
 	SHfloat max_pt[2] = {};
 	index_bounds(h_a, min_pt, max_pt);

	line_t ray(a, end);
	vec_t bbox_min(min_pt[0], min_pt[1]);
	vec_t bbox_max(max_pt[0], max_pt[1]);
	vec_t intercept;
	vec_t orig = a;

	ulong out_code     = calc_outcode(bbox_min, bbox_max, end);
	cross_type_t cross = SKEW_NO_CROSS;
	
	const double m = 0.0; // kErrorMargin; // 0.1 -> to consider larger boxes...but WHY???
	
	if (out_code == 0) {    // if end in same box as start of ray
		intercept = end;
		goto finish;
	} else {
		if (out_code & CLIP_LEFT) {
            // PPK: consider slightly larger box (m)
			line_t left(bbox_min, vec_t(bbox_min.x-m, bbox_max.y+m));
			cross = segment_intersect(left, ray, &intercept);
		} else if (out_code & CLIP_RIGHT) {
            // PPK: again consider slightly larger box (m)
			line_t right(vec_t(bbox_max.x+m, bbox_min.y-m), bbox_max);
			cross = segment_intersect(right, ray, &intercept);
		}

        // PPK: it might have already found the crossing 
		if (cross == SKEW_CROSS)
			goto finish;

		if (out_code & CLIP_TOP) {
			line_t top(bbox_min, vec_t(bbox_max.x+m, bbox_min.y-m));
			cross = segment_intersect(top, ray, &intercept);
		} else if (out_code & CLIP_BOTTOM) {
			line_t bottom(vec_t(bbox_min.x-m, bbox_max.y+m), bbox_max);
			cross = segment_intersect(bottom, ray, &intercept);
		}

		goto finish;
	}

finish:
	vec_t d(a, intercept);  // PPK: ray = start till interception (can be end of ray)
	vec_t vray(a, end);
    double lengthAEnd = vray.magnitude();   // PPK: to ensure no overshoot over end
	double l = d.magnitude() + kTraversalOffset;
	d.normalize();
    if (l > lengthAEnd){
        // std::ofstream output_file("./travRay.dat", std::ios_base::app); //HELPER
        // output_file << a.x << " " << a.y << " " 
        //             << intercept.x << " " << intercept.y << " " 
        //             << end.x << " " << end.y << " " 
        //             << min_pt[0] << " " << min_pt[1] << " " 
        //             << max_pt[0] << " " << max_pt[1] << " " 
        //             << l << " " << lengthAEnd << std::endl;
        // std::cout << "l: " <<  l << " lengthAEnd:" << lengthAEnd << std::endl;
        l = lengthAEnd;     // PPK: to ensure no overshoot over end
    }
	d.scale(l);  // PPK: thus d = (|ray|+0.0001) ray/|ray|
	a.add(d);   // PPK: till interception + 0.0001 in same direction
 
    // PPK: returns distance between start point and box-interception or
    //      end point 
	return vec_t::distance(orig, a);
}


/*PPK:
 * saves all hash-bin-ids in list "results" which are between 
 * point a and b
 */ 
void space::bins_between(SHfloat a[2], SHfloat b[2], SHuint h_b,
	std::set<SHuint> *result)
{
    // for debugging: redefine bins_between(..., SHuint c) and set in ray_casting.h to bins_between(..., SHuint c=0)
    // if (c > 500){
    //     std::cout<< "start:" << a[0] << " " << a[1] << std::endl;
    //     std::cout<< "end:" << b[0] << " " << b[1] << std::endl;
 	//     SHfloat min_pt[2],
 	//             max_pt[2];
 	//     index_bounds(h_b, min_pt, max_pt);
    //     std::cout<< "minpt:" << min_pt[0] << " " << min_pt[1] << std::endl;
    //     std::cout<< "maxpt:" << max_pt[0] << " " << max_pt[1] << std::endl;
    // }
    // c++;
	SHuint h_a = hash(a);
	if (h_a == SH_OUT_OF_BOUNDS) return;
	result->insert(h_a);
	if (h_a == h_b) return;
	
	vec_t va(a[0], a[1]),
	      vb(b[0], b[1]);
    // PPK: computes the distance from va to interception of box h_a
    //      also sets va to this interception point which is
    //      in the next box (not h_a anymore)
	double moved_by = traverse_ray(va, vb, h_a);
	a[0] = va.x;
	a[1] = va.y;
	
    // PPK: if a and b are not same point (only case with moved_by==0)
    //          --> not necessary since than h_a == h_b --> return
	if (moved_by > 0.0)
		bins_between(a, b, h_b, result);
}

/*PPK:
 * as cast_ray but modified such that it does not return hit
 * if the "ray" crosses with a line of body with "body_id"
 */
ray_result_t space::cast_ray_from_body(line_t ray, uint body_id)
{
	vec_t &start = ray.start;
	vec_t &end   = ray.end;
	
	ray_result_t result;
	result.hit_point = vec_t::kZero;
	result.hit       = false;
	result.line      = NULL;

	SHfloat s[] = { start.x, start.y };
	SHfloat e[] = { end.x,     end.y };

	SHfloat moved_by = 0.0;
	SHuint h_s;

	// keep track of bins traversed
	std::set<SHuint> traversed;

	do {
		h_s = hash(s);

		// check for any strange traversal patterns
		if (h_s == SH_OUT_OF_BOUNDS) break;
		else if (traversed.count(h_s) > 0) break;
		else traversed.insert(h_s);

		std::set<line_t*> *lines = index(h_s);
		
		if (!lines->empty()) {
			double nearest = DBL_MAX;
			for (std::set<line_t*>::iterator i = lines->begin();
			     i != lines->end(); i++)
			{
				vec_t hit_point = vec_t::kZero;
				line_t   r_copy = ray;
				line_t   c_copy = **i;
                if (c_copy.body->owner_tid != body_id)
                {
				    cross_type_t ct =
				    	segment_intersect(r_copy, c_copy, &hit_point);

				    // check for an intersection
				    if (ct == SKEW_CROSS) {
				    	SHfloat min_pt[2] = {};
				    	SHfloat max_pt[2] = {};
				    	index_bounds(h_s, min_pt, max_pt);
				    	vec_t min_v(min_pt[0], min_pt[1]),
				    	      max_v(max_pt[0], max_pt[1]);

				    	// check that the intersection is in the box
				    	if (calc_outcode(min_v, max_v, hit_point, kErrorMargin)
				    	    == 0)
				    	{
				    		double dist = vec_t::distance(ray.start, hit_point);
				    		if (dist < nearest) {
				    			nearest          = dist;
				    			result.hit_point = hit_point;
				    			result.hit       = true;
				    			result.line      = *i;

				    			point_pos_type_t pt 
				    				= classify_point(start, c_copy);

                                // PPK: compute normal vector to colliding line 
				    			if (pt == LEFT) {
				    				result.hit_normal = vec_t(
				    					c_copy.start.y - c_copy.end.y,
				    					c_copy.end.x   - c_copy.start.x
				    				);
				    			} else if (pt == RIGHT) {
				    				result.hit_normal = vec_t(
				    					c_copy.end.y   - c_copy.start.y,
				    					c_copy.start.x - c_copy.end.x
				    				);
				    			}

				    			result.hit_normal.normalize();
				    		}
				    	}
				    }
                }
			}
		}

		if (result.hit != true) {
			vec_t current(s[0], s[1]);
			vec_t end(e[0], e[1]);

            // PPK: moved_by==0 ONLY if "current" and "end" are same point 
			moved_by = traverse_ray(current, end, h_s);
			s[0] = current.x;
			s[1] = current.y;
		}
	} while (moved_by > 0.0 && result.hit != true);

	return result;
}

/*PPK:
 * as cast_ray_from_body but neglecting lines of multiple bodies
 */
ray_result_t space::cast_ray_neglect_bodies(line_t ray,
                                            std::vector<uint> &bodies_id)
{
	vec_t &start = ray.start;
	vec_t &end   = ray.end;
	
	ray_result_t result;
	result.hit_point = vec_t::kZero;
	result.hit       = false;
	result.line      = NULL;

	SHfloat s[] = { start.x, start.y };
	SHfloat e[] = { end.x,     end.y };

	SHfloat moved_by = 0.0;
	SHuint h_s;

	// keep track of bins traversed
	std::set<SHuint> traversed;

	do {
		h_s = hash(s);

		// check for any strange traversal patterns
		if (h_s == SH_OUT_OF_BOUNDS) break;
		else if (traversed.count(h_s) > 0) break;
		else traversed.insert(h_s);

		std::set<line_t*> *lines = index(h_s);
		
		if (!lines->empty()) {
			double nearest = DBL_MAX;
			for (std::set<line_t*>::iterator i = lines->begin();
			     i != lines->end(); i++)
			{
				vec_t hit_point = vec_t::kZero;
				line_t   r_copy = ray;
				line_t   c_copy = **i;
                // PPK: only checks if line is not in bodies_id  
                if (std::find(bodies_id.begin(), bodies_id.end(), 
                              c_copy.body->owner_tid) == bodies_id.end())
                {
				    cross_type_t ct =
				    	segment_intersect(r_copy, c_copy, &hit_point);

				    // check for an intersection
				    if (ct == SKEW_CROSS) {
				    	SHfloat min_pt[2] = {};
				    	SHfloat max_pt[2] = {};
				    	index_bounds(h_s, min_pt, max_pt);
				    	vec_t min_v(min_pt[0], min_pt[1]),
				    	      max_v(max_pt[0], max_pt[1]);

				    	// check that the intersection is in the box
				    	if (calc_outcode(min_v, max_v, hit_point, kErrorMargin)
				    	    == 0)
				    	{
				    		double dist = vec_t::distance(ray.start, hit_point);
				    		if (dist < nearest) {
				    			nearest          = dist;
				    			result.hit_point = hit_point;
				    			result.hit       = true;
				    			result.line      = *i;

				    			point_pos_type_t pt 
				    				= classify_point(start, c_copy);

                                // PPK: compute normal vector to colliding line 
				    			if (pt == LEFT) {
				    				result.hit_normal = vec_t(
				    					c_copy.start.y - c_copy.end.y,
				    					c_copy.end.x   - c_copy.start.x
				    				);
				    			} else if (pt == RIGHT) {
				    				result.hit_normal = vec_t(
				    					c_copy.end.y   - c_copy.start.y,
				    					c_copy.start.x - c_copy.end.x
				    				);
				    			}

				    			result.hit_normal.normalize();
				    		}
				    	}
				    }
                }
			}
		}

		if (result.hit != true) {
			vec_t current(s[0], s[1]);
			vec_t end(e[0], e[1]);

            // PPK: moved_by==0 ONLY if "current" and "end" are same point 
			moved_by = traverse_ray(current, end, h_s);
			s[0] = current.x;
			s[1] = current.y;
		}
	} while (moved_by > 0.0 && result.hit != true);

	return result;
}


}; // rc
