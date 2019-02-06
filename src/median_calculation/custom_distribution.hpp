/*
This file is part of UMAP.  For copyright information see the COPYRIGHT
file in the top level directory, or at
https://github.com/LLNL/umap/blob/master/COPYRIGHT
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License (as published by the Free
Software Foundation) version 2.1 dated February 1999.  This program is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the terms and conditions of the GNU Lesser General Public License
for more details.  You should have received a copy of the GNU Lesser General
Public License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef UMAP_APPS_MEDIAN_CALCULATION_CUSTOM_DISTRIBUTION_HPP
#define UMAP_APPS_MEDIAN_CALCULATION_CUSTOM_DISTRIBUTION_HPP

#include <random>

namespace median {
	
	//function for sampling a random value from the inverse cdf (inversion is done implicitly)
	float cdf_sample(float *cdf, const uint32_t &nbins, const float &minBound, const float maxBound)
	{
		float r = drand48();
		float *ptr = std::lower_bound(cdf, cdf + nbins + 1, r);
		int off = std:max(0, (int)(ptr - cdf - 1));
		float t = (r -cdf[off])/(cdf[off + 1] - cdf[off]);
		float x = (off + t)/(float)(nbins);
		return minBound + (maxBound-minBound)*x;
	}

	//somehow load in pdf
	//get min_bound, max_bound, and num_bins as values from pdf_bins
	//get pdf values and put into array
	
	int nbins = ???;
	float minBound = ???, maxBound = ???;
	float pdf[nbins] = ???;

	//create cdf from pdf
	float cdf[nbins+1], dx = (maxBound - minBound)/nbins, sum = 0;
	cdf[0] = 0.f;
	for (int n= 1; n<nbins; ++n) 
	{
		cdf[n] = cdf[n-1] + pdf[n-1];
		sum += pdf[n-1];
	}
	cdf[nbins] = 1;


} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_CUSTOM_DISTRIBUTION_HPP
