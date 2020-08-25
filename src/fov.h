/*
 *  fov.h
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

#ifndef FOV_RAY_SPACE_H
#define FOV_RAY_SPACE_H

#include <vector>
// #include <opencv/cv.h>
#include "ray_casting.h"

#ifndef uint
typedef unsigned int uint;
#endif

namespace fov {

/*
 *	Creates an efficient data structure for storing geometry and querying the
 *	space for ray intersections.
 */
class ray_space {
public:
	ray_space(gm::vec_t dim);
	~ray_space();

	rc::vec_t   dim();
	rc::space*  space();
	rc::body_t* boundary();
	void        reset();

	// creating and adding geometry data
	rc::body_t* create_body_for_boundary(uint owner_tid);
	rc::body_t* create_body_for_rectangle(
		rc::vec_t x0, rc::vec_t x1, rc::vec_t x2, rc::vec_t x3,
		uint owner_tid
	);
	rc::body_t* create_body_for_circle(
		rc::vec_t center, double r, uint f, uint owner_tid=-1
	);
	std::vector<rc::body_t*> create_body_for_circular_detector(
		rc::vec_t center, double r, uint f, uint &tid_iterator
	);
	// rc::body_t* create_body_from_sequence(CvSeq* seq, uint owner_tid=-1);
	rc::body_t* create_body_from_sequence(
		std::vector<rc::line_t> &lines, uint owner_tid=-1
	);
	void add_body_to_space(rc::body_t *body);

	// querying stored geometry for ray intersections
	rc::vec_t terminus_of_ray(
		rc::vec_t orig, rc::vec_t end, double project=0.0,
		rc::body_t **hitb=NULL
	);
	rc::vec_t terminus_of_ray_from_body(
		rc::vec_t orig, rc::vec_t end, uint body_id,
        double project=0.0, rc::body_t **hitb=NULL
	);
	rc::vec_t terminus_of_ray_neglect_bodies(
		rc::vec_t orig, rc::vec_t end,
        std::vector<uint> &bodies_id,
        double project=0.0, rc::body_t **hitb=NULL
	);
	void view_at_point(
		rc::vec_t orig, double theta_offset, double winding,
		rc::vec_t *termini, rc::body_t **hit_body, uint res
	);
	void view_at_point_from_body(
		rc::vec_t orig, double theta_offset, double winding,
		rc::vec_t *termini, rc::body_t **hit_body, uint res,
        uint body_id
	);
	void view_at_point_neglect_bodies(
		rc::vec_t orig, double theta_offset, double winding,
		rc::vec_t *termini, rc::body_t **hit_body, uint res,
        std::vector<uint> &bodies_id
	);
	bool can_see(
		rc::vec_t orig, rc::vec_t pt, double tolerance,
		bool exclude_boundary = false
	);

private:
	rc::space  *m_space;
	rc::vec_t  m_dim;
	rc::body_t *m_boundary;
};


}; // fov

#endif // FOV_RAY_SPACE_H
