/*  SocialForces
    defines social force models which describe the interactions
    between agents in SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#include "social_forces.h"

// Social-Force-Models------------------------------------------------
// Distance Regulation and Alignment (DRA)
void SFM_DRA(params *ptrSP, std::vector<double> &u_ji,
             double dist_interaction, std::vector<double> &v_ji,
             std::vector<double> &f0, std::vector<double> &f1){
    /* computes the 4 forces assuming input of c = {0, 0, 0, 0}
     * and fx = {0, 0}
     */
    double fstrength = 0.0;
    // FORCE0 (distance regulation = repulsion + attraction)
    fstrength = (2 * SigmoidFunction(dist_interaction, ptrSP->rep_range,
                                ptrSP->rep_steepness)) - 1;
    f0 = vec_mul(u_ji, -fstrength);
    // FORCE1 (alignment)
    f1 = v_ji;
}
