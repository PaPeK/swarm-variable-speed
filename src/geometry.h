/*
 *  geometry.h
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

#ifndef FOV_GEOMETRY_H
#define FOV_GEOMETRY_H

#include <stdlib.h>

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif

namespace gm {


class vec_t {
public:
	double x;
	double y;

	vec_t();
	vec_t(double x_, double y_);
	vec_t(double *x_);
	vec_t(vec_t &u, vec_t &v);
	vec_t(double *u, double *v);

	void add(vec_t &v);
	void sub(vec_t &v);
	void scale(double s);
	void scale(vec_t &v);
	void normalize();
	void rotate(double theta);
	double magnitude();

	static double distance(vec_t &a, vec_t &b);
	static double angle_between(vec_t &a, vec_t &b);
	static double dot(vec_t &a, vec_t &b);

	vec_t scaled(double s);
	vec_t orthogonal();

	static const vec_t kZero;
	
	friend bool operator==(const vec_t &v, const vec_t &u);
	friend bool operator!=(const vec_t &v, const vec_t &u);

	friend vec_t operator+(const vec_t &v, const vec_t &u);
	friend vec_t operator-(const vec_t &v, const vec_t &u);
	friend vec_t operator*(const vec_t &v, const vec_t &u);
	friend vec_t operator*(const vec_t &v, const double m);
	friend vec_t operator/(const vec_t &v, const vec_t &u);
	friend vec_t operator/(const vec_t &v, const double m);
};


class line_t;

class body_t {
public:
	uint bid;
	uint owner_tid;
	uint total_lines;
	line_t *lines;      // PPK: is that save? doesn't it cause seg-faults -> use vector or new ??? ->>>> In combination with space::add_body it works fine (use of new)

	body_t();

	void add_line(line_t line);
	
private:
	uint m_next_line;
};


class line_t {
public:
	vec_t  start;
	vec_t  end;
	body_t *body;

	line_t()
		: start(), end(), body(NULL)
	{ }

	line_t(vec_t s, vec_t e, body_t *b=NULL)
		: start(s), end(e), body(b)
	{ }

	// vec_t  s(double t);
	line_t scaled(double scale);
	line_t extended(double amount);
};


class transform_t {
public:
	transform_t();
	transform_t(vec_t &axis);
	transform_t(vec_t &a, vec_t &b);
	transform_t(vec_t &offset, vec_t &axis, double &scale);
	transform_t(vec_t offset, double scale, double theta);

	void apply(vec_t &v);
	void unapply(vec_t &v);

	gm::vec_t& offset();
	double scale();
	double theta();

private:
	vec_t  m_offset;
	double m_scale;
	double m_theta;
};


}; // gm

#endif // FOV_GEOMETRY_H
