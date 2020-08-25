/*
 *  fov.cpp
 *  This file is part of FOVEA
 *
 *  Created by Colin Twomey on 11.18.2010
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

#include "fov.h"

#define DEG_TO_RAD(x) (x * M_PI / 180.0)
#define RAD_TO_DEG(x) (x * 180.0 / M_PI)
#ifndef SQ
#define SQ(a)         ((a)*(a))
#endif

namespace fov {

/*
 *	Setup a spatial-hash data structure for efficient ray-to-object collision
 *	detection.
 */
ray_space::ray_space(gm::vec_t dim)
{
	// define the structure of the spatial hash
	m_dim = dim;

	// create and initalize the spatial hash
	SHuint  bdim[]   = { 16, 16 };
	SHfloat size[]   = { dim.x, dim.y };
	SHfloat offset[] = { 0, 0 };

	m_space = new rc::space(bdim, size, offset);
	m_boundary = NULL;
}


ray_space::~ray_space(void)
{
	m_space->clear_bodies();
	delete m_space;
}


rc::vec_t ray_space::dim(void)
{
	return m_dim;
}


rc::space* ray_space::space(void)
{
	return m_space;
}


rc::body_t* ray_space::boundary(void)
{
	return m_boundary;
}


void ray_space::reset(void)
{
	m_space->clear_bodies();
	m_space->clear();
}


rc::body_t* ray_space::create_body_for_boundary(uint owner_tid)
{
	m_boundary = m_space->add_body(4, owner_tid);
	rc::vec_t origin   = rc::vec_t(       0.0,        0.0);
	rc::vec_t opposite = rc::vec_t(   m_dim.x,    m_dim.y);
	rc::vec_t opx_ory  = rc::vec_t(opposite.x,   origin.y);
	rc::vec_t orx_opy  = rc::vec_t(  origin.x, opposite.y);
	
	m_boundary->add_line(rc::line_t( origin,  opx_ory));
	m_boundary->add_line(rc::line_t(orx_opy, opposite));
	
	m_boundary->add_line(rc::line_t( origin,  orx_opy));
	m_boundary->add_line(rc::line_t(opx_ory, opposite));

	return m_boundary;
}


rc::body_t* ray_space::create_body_for_rectangle(
	rc::vec_t x0, rc::vec_t x1, rc::vec_t x2, rc::vec_t x3,
	uint owner_tid)
{
	rc::body_t *body = m_space->add_body(4, owner_tid);

    body->add_line(rc::line_t(x0, x1));
	body->add_line(rc::line_t(x1, x2));
	body->add_line(rc::line_t(x2, x3));
	body->add_line(rc::line_t(x3, x0));

	return body;
}


rc::body_t* ray_space::create_body_for_circle(
	rc::vec_t center, double r, uint f, uint owner_tid)
{
	rc::body_t *body = m_space->add_body(f, owner_tid);

	const double kDtheta = (2.0 * M_PI) / f;
	rc::vec_t start  = rc::vec_t(center.x, center.y + r);
	rc::vec_t prev_p = start;

	for (uint i = 1; i < f; i++) {
		double theta = i * kDtheta;
		rc::vec_t p(sin(theta), cos(theta));
		p.scale(r);
		p.add(center);

		body->add_line(rc::line_t(prev_p, p));
		prev_p = p;
	}

	// close the polygon
	body->add_line(rc::line_t(prev_p, start));

	return body;
}


std::vector<rc::body_t*> ray_space::create_body_for_circular_detector(
	rc::vec_t center, double r, uint f, uint &tid_iterator)
{
	std::vector<rc::body_t*> bodies(f);

	const double kDtheta = (2.0 * M_PI) / f;
	rc::vec_t start  = rc::vec_t(center.x, center.y + r);
	rc::vec_t prev_p = start;

	for (uint i = 1; i < f; i++) {
		rc::body_t *body = m_space->add_body(1, tid_iterator++);

		double theta = i * kDtheta;
		rc::vec_t p(sin(theta), cos(theta));
		p.scale(r);
		p.add(center);

		body->add_line(rc::line_t(prev_p, p));
		prev_p = p;

		bodies[i-1] = body;
	}

	// close the polygon
	rc::body_t *body = m_space->add_body(1, tid_iterator++);
	body->add_line(rc::line_t(prev_p, start));
	bodies[f-1] = body;

	return bodies;
}


// /*
//  *	Given an OpenCV line sequence, generate a set of static line segment
//  *	shapes and attach them to a body representing the whole sequence. Note
//  *  that the body is destroyed by the ray casting subsystem, and should
//  *  never be released directly by the caller.
//  */
// rc::body_t* ray_space::create_body_from_sequence(CvSeq* seq, uint owner_tid)
// {
// 	rc::body_t *body = m_space->add_body(seq->total, owner_tid);
// 
// 	for (int i = 0; i < seq->total; i++) {
// 		CvPoint *p = (CvPoint*)cvGetSeqElem(seq, i);
// 		CvPoint *q = (CvPoint*)cvGetSeqElem(seq, (i+1) % seq->total);
// 
// 		rc::vec_t  v0 = rc::vec_t(p->x, p->y);
// 		rc::vec_t  v1 = rc::vec_t(q->x, q->y);
// 
// 		body->add_line(rc::line_t(v0, v1));
// 	}
// 
// 	return body;
// }

/*
 *  Creates a body from an array of line segments. Note that the body is
 *  destroyed by the ray casting subsystem, and should never be released
 *  directly by the caller.
 */
rc::body_t* ray_space::create_body_from_sequence(
	std::vector<rc::line_t> &lines,
	uint owner_tid)
{
	rc::body_t *body = m_space->add_body(lines.size(), owner_tid);

	for (uint i = 0; i < lines.size(); i++) {
		rc::line_t l = lines[i];
		body->add_line(l);
	}

	return body;
}


/*
 *	Insert the given body and its shapes into the spatial hash.
 */
void ray_space::add_body_to_space(rc::body_t *body)
{
	m_space->attach_body(body);
}

/*
 *	Cast a ray from orig to end, and return the point at which it collides
 *	with any shape in the spatial hash (or end, if it doesn't intersect
 *	anything).  An offset in the direction of the normal of the surface hit
 *	can be specified using the project parameter.
 */
rc::vec_t ray_space::terminus_of_ray(
	rc::vec_t  orig,
	rc::vec_t  end,
	double     project,
	rc::body_t **hitb)
{
	rc::ray_result_t rr = m_space->cast_ray(rc::line_t(orig, end));
	if (rr.hit) {
		if (hitb) *hitb = rr.line->body;
		rc::vec_t hit = rr.hit_point;
		rc::vec_t t   = rc::vec_t(hit.x + rr.hit_normal.x * project,
		                          hit.y + rr.hit_normal.y * project);
		return t;
	} else if (hitb) *hitb = NULL;

	return end;
}

/*
 * PPK: same as terminus_of_ray but the ray ignores lines of 
 *      body with "body_id"
 */
rc::vec_t ray_space::terminus_of_ray_from_body(
	rc::vec_t  orig,
	rc::vec_t  end,
    uint body_id,
	double     project,
	rc::body_t **hitb)
{
	rc::ray_result_t rr = m_space->cast_ray_from_body(rc::line_t(orig, end),
                                                      body_id);
	if (rr.hit) {
		if (hitb) *hitb = rr.line->body;
		rc::vec_t hit = rr.hit_point;
		rc::vec_t t   = rc::vec_t(hit.x + rr.hit_normal.x * project,
		                          hit.y + rr.hit_normal.y * project);
		return t;
	} else if (hitb) *hitb = NULL;

	return end;
}

/*
 * PPK: same as terminus_of_ray but the ray ignores lines of 
 *      body with "body_id"
 */
rc::vec_t ray_space::terminus_of_ray_neglect_bodies(
	rc::vec_t  orig,
	rc::vec_t  end,
    std::vector<uint> &bodies_id,
	double     project,
	rc::body_t **hitb)
{
	rc::ray_result_t rr = m_space->cast_ray_neglect_bodies(rc::line_t(orig, end),
                                                           bodies_id);
	if (rr.hit) {
		if (hitb) *hitb = rr.line->body;
		rc::vec_t hit = rr.hit_point;
		rc::vec_t t   = rc::vec_t(hit.x + rr.hit_normal.x * project,
		                          hit.y + rr.hit_normal.y * project);
		return t;
	} else if (hitb) *hitb = NULL;

	return end;
}

/*
 *	Calculates a 360° view at a given point and with a given resolution
 *	(where the resolution just determines the number of rays to test).
 */
void ray_space::view_at_point(
	rc::vec_t  orig,
	double     theta_offset,
	double     winding,
	rc::vec_t  *termini,
	rc::body_t **hit_body,
	uint       res)
{
	const double kMaxLength = sqrt(SQ(m_dim.x) + SQ(m_dim.y));

	double d_theta = (1.0 / (double)res) * M_PI * 2.0; // 2pi for full, pi for half!
	// NOTE: at one point this was set to just pi, but not clear why this was the case.
	d_theta = copysign(d_theta, winding);

	for (uint i = 0; i < res; i++) {
		double theta = i * d_theta + theta_offset;
		rc::vec_t end = rc::vec_t(orig.x + sin(theta) * kMaxLength,
		                          orig.y + cos(theta) * kMaxLength);
		termini[i] = terminus_of_ray(orig, end, 0.0, &hit_body[i]);
	}
}

/*
 * PPK: same as view_at_point but ignores lines which belong to
 *      body with "body_id"
 */
void ray_space::view_at_point_from_body(
	rc::vec_t  orig,
	double     theta_offset,
	double     winding,
	rc::vec_t  *termini,
	rc::body_t **hit_body,
	uint       res,
    uint       body_id)
{
	const double kMaxLength = sqrt(SQ(m_dim.x) + SQ(m_dim.y));

	double d_theta = (1.0 / (double)res) * M_PI * 2.0; // 2pi for full, pi for half!
	// NOTE: at one point this was set to just pi, but not clear why this was the case.
	d_theta = copysign(d_theta, winding);

	for (uint i = 0; i < res; i++) {
		double theta = i * d_theta + theta_offset;
		rc::vec_t end = rc::vec_t(orig.x + sin(theta) * kMaxLength,
		                          orig.y + cos(theta) * kMaxLength);
		termini[i] = terminus_of_ray_from_body(orig, end, body_id, 0.0, 
                                               &hit_body[i]);
	}
}

/*
 * PPK: same as view_at_point_from_bofy but ignores lines 
 *      from multiple bodies which are in "bodies_id"
 */
void ray_space::view_at_point_neglect_bodies(
	rc::vec_t  orig,
	double     theta_offset,
	double     winding,
	rc::vec_t  *termini,
	rc::body_t **hit_body,
	uint       res,
    std::vector<uint> &bodies_id)
{
	const double kMaxLength = sqrt(SQ(m_dim.x) + SQ(m_dim.y));

	double d_theta = (1.0 / (double)res) * M_PI * 2.0; // 2pi for full, pi for half!
	// NOTE: at one point this was set to just pi, but not clear why this was the case.
	d_theta = copysign(d_theta, winding);

	for (uint i = 0; i < res; i++) {
		double theta = i * d_theta + theta_offset;
		rc::vec_t end = rc::vec_t(orig.x + sin(theta) * kMaxLength,
		                          orig.y + cos(theta) * kMaxLength);
		termini[i] = terminus_of_ray_neglect_bodies(orig, end, bodies_id, 0.0, 
                                                    &hit_body[i]);
	}
}

/*
 *  Determine if a given point (pt) is visible from a ray started at the given
 *  origin (orig).  Return true if yes, false otherwise
 */
bool ray_space::can_see(
	rc::vec_t orig,
	rc::vec_t pt,
	double tolerance,
	bool exclude_boundary)
{
	rc::body_t *hit    = NULL;
	rc::vec_t rays_end = terminus_of_ray(orig, pt, 0.0, &hit);

	return (hit == NULL) || (exclude_boundary && hit == m_boundary);
}


}; // fov
