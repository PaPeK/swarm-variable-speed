/*
 *  geometry.cpp
 *  This file is part of FOVEA
 *
 *  Created by Colin Twomey on 12.23.2011
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

#include <cmath>
#include "geometry.h"

#ifndef SQ
#define SQ(a) ((a)*(a))
#endif

namespace gm {


const vec_t vec_t::kZero(0.0, 0.0);


vec_t::vec_t(void)
	: x(0.0), y(0.0)
{ }


vec_t::vec_t(double x_, double y_)
	: x(x_), y(y_)
{ }


vec_t::vec_t(double *x_)
	: x(x_[0]), y(x_[1])
{ }


vec_t::vec_t(vec_t &u, vec_t &v)
{
	x = v.x - u.x;
	y = v.y - u.y;
}


vec_t::vec_t(double *u, double *v)
{
	x = v[0] - u[0];
	y = v[1] - u[1];
}


void vec_t::add(vec_t &v)
{
	x += v.x;
	y += v.y;
}


void vec_t::sub(vec_t &v)
{
	x -= v.x;
	y -= v.y;
}


void vec_t::scale(double s)
{
	x *= s;
	y *= s;
}


void vec_t::scale(vec_t &v)
{
	x *= v.x;
	y *= v.y;
}


void vec_t::normalize(void)
{
	double m = 1.0f / magnitude();
	x *= m;
	y *= m;
}


double vec_t::magnitude(void)
{
	return sqrt(SQ(x) + SQ(y));
}

/*
 *	Rotates the 2d vector by a specified angle theta.
 */
void vec_t::rotate(double theta)
{
	vec_t d;

	double sin_theta = sinf(theta),
	       cos_theta = cosf(theta);

	d.x = x * cos_theta - y * sin_theta;
	d.y = x * sin_theta + y * cos_theta;

	x = d.x;
	y = d.y;
}


double vec_t::distance(vec_t &a, vec_t &b)
{
	return sqrt(SQ(b.x - a.x) + SQ(b.y - a.y));
}

/*
 *	Returns the signed angle on the interval [-pi, pi] between vector a and b.
 *  WARNING: assumes vectors a and b are normal!
 */
double vec_t::angle_between(vec_t &a, vec_t &b)
{
	double d = a.x * b.x + a.y * b.y;
	if		(d < -1.0f) d = -1.0f;
	else if	 (d > 1.0f) d =  1.0f;
	return copysign(acos(d), a.y * b.x - a.x * b.y);
}


double vec_t::dot(vec_t &a, vec_t &b)
{
	return a.x * b.x + a.y * b.y;
}


vec_t vec_t::scaled(double s)
{
	return vec_t(x*s, y*s);
}


vec_t vec_t::orthogonal(void)
{
	return vec_t(-y, x);
}


bool operator==(const vec_t &v, const vec_t &u)
{
	return (v.x == u.x) && (v.y == u.y);
}


bool operator!=(const vec_t &v, const vec_t &u)
{
	return !(v == u);
}


vec_t operator+(const vec_t &v, const vec_t &u)
{
	return vec_t(v.x + u.x, v.y + u.y);
}


vec_t operator-(const vec_t &v, const vec_t &u)
{
	return vec_t(v.x - u.x, v.y - u.y);
}


vec_t operator*(const vec_t &v, const vec_t &u)
{
	return vec_t(v.x * u.x, v.y * u.y);
}


vec_t operator*(const vec_t &v, const double m)
{
	return vec_t(v.x * m, v.y * m);
}


vec_t operator/(const vec_t &v, const vec_t &u)
{
	return vec_t(v.x / u.x, v.y / u.y);
}


vec_t operator/(const vec_t &v, const double m)
{
	return vec_t(v.x / m, v.y / m);
}


// #pragma mark -

body_t::body_t()
	: bid(-1),
	  owner_tid(-1),
	  total_lines(0),
	  m_next_line(0),
	  lines(NULL)
{ }


void body_t::add_line(line_t line)
{
	if (lines == NULL) return;
	this->lines[m_next_line]      = line;
	this->lines[m_next_line].body = this;
	this->m_next_line++;
}


// #pragma mark -

/*
 *  Returns the coordinate along the line length t * distance(start, end) from
 *  the start.  Eg. at t = 0, return start, and at t = 1 returns end.
 */
// vec_t line_t::s(double t)
// {
// 	return gm::vec_t(start.x * t + end.x * (1.0-t),
// 	                 start.y * t + end.y * (1.0-t))
// }


line_t line_t::scaled(double scale)
{
	return line_t(start.scaled(scale), end.scaled(scale));
}


line_t line_t::extended(double amount)
{
	vec_t dir(start, end);
	dir.normalize();
	dir.scale(amount);

	vec_t s = start;
	vec_t e = end;

	e.add(dir);
	dir.scale(-1);
	s.add(dir);

	return line_t(s, e);
}


// #pragma mark -

transform_t::transform_t(void)
{
	m_theta  = 0.0;
	m_scale  = 1.0;
	m_offset = gm::vec_t(0.0, 0.0);
}


transform_t::transform_t(vec_t &axis)
{
	static vec_t kXAxis(1, 0);
	axis.normalize();

	m_theta  = vec_t::angle_between(axis, kXAxis);
	m_scale  = 1.0;
	m_offset = gm::vec_t(0.0, 0.0);
}


transform_t::transform_t(vec_t &a, vec_t &b)
{
	static vec_t kXAxis(1, 0);

	vec_t axis(b.x - a.x,
	           b.y - a.y);
	axis.normalize();

	m_theta  = vec_t::angle_between(axis, kXAxis);
	m_scale  = 1.0;
	m_offset = a;
}


transform_t::transform_t(vec_t &offset, vec_t &axis, double &scale)
{
	static vec_t kXAxis(1, 0);
	axis.normalize();

	m_theta  = vec_t::angle_between(axis, kXAxis);
	m_scale  = scale;
	m_offset = offset;
}


transform_t::transform_t(vec_t offset, double scale, double theta)
{
	m_offset = offset;
	m_scale  = scale;
	m_theta  = theta;
}


void transform_t::apply(vec_t &v)
{
	// rotate
	v.rotate(m_theta);

	// scale and translate
	v.x = v.x * m_scale + m_offset.x;
	v.y = v.y * m_scale + m_offset.y;
}


void transform_t::unapply(vec_t &v)
{
	// undo translation and scale
	v.x = (v.x - m_offset.x) / m_scale;
	v.y = (v.y - m_offset.y) / m_scale;

	// undo rotation
	v.rotate(-m_theta);
}


gm::vec_t& transform_t::offset(void)
{
	return m_offset;
}


double transform_t::scale(void)
{
	return m_scale;
}


double transform_t::theta(void)
{
	return m_theta;
}


}; // gm
