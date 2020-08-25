/*
 *  spatial_hash.h
 *  This file is part of FOVEA
 *
 *  Created by Colin Twomey on 05.16.2010
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

#ifndef FOV_SPATIAL_HASH_H
#define FOV_SPATIAL_HASH_H

#include <float.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <set>
#include <list>
#include <utility>

#define SH_OUT_OF_BOUNDS 0xFFFF


#ifndef SHuint
typedef unsigned int SHuint;
#endif

#ifndef SHull
typedef unsigned long long SHull;
#endif

#ifndef SHfloat
typedef double SHfloat;
#endif


template<SHuint D, typename T>
class spatial_hash {
public:
	/*
	 *	Construct a spatial hash with the given number of dimensions (D), the
	 *	number of bins per dimension (dim[D]), and the mapping of real space
	 *	to the bins (size[D] space per bin). 
     *	PPK: I think size[D] is the size of space considered and m_size is space per  
     *	     bin (which is computed in the code)
	 */
	spatial_hash(SHuint dim[D], SHfloat size[D], SHfloat offset[D])
	{
		m_update = 0ULL;
		
		memcpy(m_dim, dim, sizeof(T) * D);  // PPK: why sizeof(T)... it should be sizeof(SHuint) AND why copy dim to protected variable m_dim.... because other memberfctn. can access it
		remap(size, offset);
		
		m_big_dim = 1;      // # of total bins
//		#pragma unroll D
		for (SHuint i = 0; i < D; i++)
			m_big_dim *= m_dim[i];
		
		m_hash = new std::set<T>[m_big_dim];
		
		m_hash_update = new SHuint[m_big_dim];
		memset(m_hash_update, m_update, sizeof(SHuint) * m_big_dim);
	}

	/*
	 *	Destroy the spatial hash.
	 */
	~spatial_hash(void)
	{
		delete [] m_hash;
		delete [] m_hash_update;
	}

	/*
	 *	Preserve the spatial hash structure itself, but clean out any used
	 *	bins to make a clean slate for future use.
	 */
	void clear(void)
	{
		m_update += 1ULL;
		if (m_update == ~0ULL) {
			m_update = 0ULL;
			for (SHuint i = 0; i < m_big_dim; i++) {
				m_hash[i].clear();
				m_hash_update[i] = m_update;
			}
		}
	}

	/*
	 *	Change the real-space to bin-space mapping.
	 */
	void remap(SHfloat size[D], SHfloat offset[D])
	{
//		#pragma unroll D
		for (SHuint i = 0; i < D; i++)
			m_size[i] = size[i] / (SHfloat)m_dim[i];    // PPK: why?? m_size=sizeOfBin/#ofBin -> probably wrong described and size[D] is size along dimension D
		memcpy(m_offset, offset, sizeof(SHfloat) * D);  // PPK: why *2... it should be *D.... right?
	}


	SHuint dim(void)
	{
		return m_big_dim;
	}


	SHuint total_bins(void)
	{
     	return m_big_dim;
	}


	SHuint dim(SHuint d)
	{
		return m_dim[d];
	}


	SHfloat size(SHuint d)
	{
		return m_size[d];
	}

	/*
	 *	Returns the hash bin for the given spatial position ( + offset, used
	 *	to ensure pos[D] is non-negative).
	 */
	SHuint hash(SHfloat pos[D])
	{
     	SHuint h[D];
//      	#pragma unroll D
     	for (SHuint i = 0; i < D; i++) {
			SHfloat p = pos[i] + m_offset[i];
			if (p < 0.0 || p > (m_size[i] * m_dim[i]))
				return SH_OUT_OF_BOUNDS;
           	h[i] = floorf(p / m_size[i]);
			if (h[i] >= m_dim[i])
				return SH_OUT_OF_BOUNDS;
		}
          
        SHuint hash_i = 0;
        SHuint di = 1;
//         #pragma unroll D
        for (SHuint i = 0; i < D; i++) {
			hash_i += h[i] * di;
			di *= m_dim[i];
		}
		
		if (hash_i >= m_big_dim)
			return SH_OUT_OF_BOUNDS;
        return hash_i;
	}

	/*
	 *	Insert a given value into the hash bin at position hash_i.
	 */
	void place(SHuint hash_i, T value)
	{
		if (hash_i == SH_OUT_OF_BOUNDS) return;
		if (m_hash_update[hash_i] != m_update)
			m_hash[hash_i].clear();
		m_hash[hash_i].insert(value);
		m_hash_update[hash_i] = m_update;
	}

	/*
	 *	Return the values hashed into bin hash_i.
	 */
	std::set<T>* index(SHuint hash_i)
	{
		if (   hash_i == SH_OUT_OF_BOUNDS
		    || m_hash_update[hash_i] != m_update)
		{
			return &m_empty;
		}
		return &m_hash[hash_i];
	}


	void index_bounds(SHuint hash_i, SHfloat min_pt[D], SHfloat max_pt[D])
	{
		any_index_bounds(hash_i, min_pt, max_pt, m_dim, m_size, m_offset);
	}
	
	
protected:
	std::set<T> *m_hash,
	            m_empty;
	SHfloat	m_size[D],
	        m_offset[D];
	SHuint  m_big_dim,
	        m_dim[D],
			*m_hash_update;
	SHull   m_update;

	void any_index_bounds(SHuint hash_i, SHfloat min_pt[D], SHfloat max_pt[D],
	                      SHuint dim[D], SHfloat size[D], SHfloat offset[D])
	{
		if (hash_i == SH_OUT_OF_BOUNDS) return;
		SHuint di = 1;
// 		#pragma unroll D
		for (SHuint i = 0; i < D; i++) {
			SHuint h = (hash_i / di) % dim[i];
			di *= dim[i];

			min_pt[i] = h * size[i] - offset[i];
			max_pt[i] = min_pt[i] + size[i];
		}
	}
};


#endif // FOV_SPATIAL_HASH_H
