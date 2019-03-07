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

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>

namespace median {
	

class custom_distribution {
 public:
  custom_distribution(const std::string &pdf_filename)
  {
	//probability density function should come from file provided via DECam_vectors.py
	//should be perp_bins,perp_vals,para_bins,para_vals
	
	std::vector<double> perp_pdf, perp_bins, para_pdf, para_bins;
	int check = 0;
	  
	std::ifstream ifile(pdf_filename);
	
	if (!ifile.is_open()) {
	      std::cerr << "Cannot open " << pdf_filename << std::endl;
	      std::abort();
	}
	std::string line;
	std::getline(ifile, line); //skip first row of header info
	
	while (std::getline(ifile, line))  
	{
		std::istringstream iss{line};
		
		double perp_b, perp_v, para_b, para_v;
		char vala, valb, valc, vald;
		int id;
		do{
		if( iss >> id >> vala >> perp_b >> valb >> perp_v >> valc >> para_b >> vald >> para_v)
		{
		perp_bins.push_back(perp_b);
		perp_pdf.push_back(perp_v);
		para_bins.push_back(para_b);
		para_pdf.push_back(para_v);
		}}while(!iss.eof());
		
		++check;
	}
	
	nbins = check;
	para_minBound = para_bins.front(), para_maxBound = para_bins.back();
	perp_minBound = perp_bins.front(), perp_maxBound = perp_bins.back();
	para_cdf = gen_cdf(para_pdf,nbins);  
	perp_cdf = gen_cdf(perp_pdf,nbins);
  }	
  template <typename rnd_engine>
  std::vector<double> operator()(rnd_engine &engine) {
	double para = cdf_sample(engine,para_cdf,nbins,para_minBound,para_maxBound);
	double perp = cdf_sample(engine,perp_cdf,nbins,perp_minBound,perp_maxBound);
	std::vector<double> ra_dec = ecliptic_to_equatorial(para,perp);
	//convert from ra/dec to pixel coordinates
	//Should generalize in the future to be based on fits header wcs
	std::vector<double> xy_coords = {-ra_dec[0]/0.27,ra_dec[1]/0.27};
	return xy_coords;
  }
	
 private:

  std::vector<double> perp_cdf, para_cdf;
  int nbins;
  double para_minBound, perp_minBound, para_maxBound, perp_maxBound;
  
  // Function for sampling a random value from the inverse cdf (inversion is done implicitly)
  // This is slightly more complicated to allow for more precise sampling (between given cdf values)
  template <typename rnd_engine>
  double cdf_sample(rnd_engine &engine, std::vector<double> cdf, const uint32_t &nbins, const double &minBound, const double maxBound)
  {
  	std::uniform_real_distribution<double> dis(0.0, 1.0);
	double r = dis(engine);
	std::vector<double>::iterator lwr = std::lower_bound(cdf.begin(), cdf.end(), r);
  	int under_bound = std::max(0, (int)(lwr - 1 - cdf.begin()));
  	double t = (r -cdf[under_bound])/(cdf[under_bound+1] - cdf[under_bound]);
  	double x = (under_bound + t)/(double)(nbins);
  	return minBound + (maxBound-minBound)*x;
  }
  
  //function for creating a cumulative distribution function from a probability distribution function
  std::vector<double> gen_cdf(std::vector<double> pdf, int nbins)
  {
	std::vector<double> cdf(nbins);
	cdf[0] = 0.0;
	for (int n= 1; n<nbins; ++n) 
	{
		cdf[n] = cdf[n-1] + pdf[n-1];
	}
	cdf[nbins] = 1.0;
	return cdf;
  }

  //function for converting ecliptic coordinates to equatorial
  //takes in arcsecs/sec and outputs arcsecs/sec
  //assumes same epoch for both frames
  std::vector<double> ecliptic_to_equatorial(double para, double perp)
  {
	double ecl_obl = 0.4090926; //radians
	double dec = asin(sin(ecl_obl)*sin(para/206265.0)*cos(perp/206265.0)+cos(ecl_obl)*sin(perp/206265.0))*206265.0;
	double ra = atan2((cos(ecl_obl)*sin(para/206265.0) - sin(ecl_obl)*tan(perp/206265.0)),cos(para/206265.0))*206265.0;
	std::vector<double> v = {ra,dec}; //now in arcsec/sec
	return v; 
  }
};	
	
} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_CUSTOM_DISTRIBUTION_HPP
