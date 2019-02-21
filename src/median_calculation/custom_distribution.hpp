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


namespace median {
	

class custom_distribution {
 public:
  custom_distribution(const std::string &pdf_filename)
  {
	//probability density function should come from file provided via DECam_vectors.py
	//should be perp_bins,perp_vals,para_bins,para_vals
	
	std::vector<double> perp_pdf, perp_bins, para_pdf, para_bins;
	int check = 0;
	  
	std::ifstream ifile(pdf_file_name);
	if (!ifile.is_open()) {
	      std::cerr << "Cannot open " << pdf_file_name << std::endl;
	      std::abort();
	}
	std::string line;
	while (std::getline(ifile, line))  
	{
		if (check==0) continue; //for the first row having column names
		++check;
		std::istrinstream iss{line};
		

		std::string temp;
		double perp_b, perp_v, para_b, para_v;
		std::getline(iss, temp, ',') >> perp_b;
		std::getline(iss, temp, ',') >> perp_v;
		std::getline(iss, temp, ',') >> para_b;
		std::getline(iss, temp) >> para_v;
		perp_bins.push_back(perp_b);
		perp_pdf.push_back(perp_v);
		para_bins.push_back(para_b);
		para_pdf.push_back(para_v);
	}

	nbins = check-1;
	para_minBound = para_bins.front(), para_maxBound = para_bins.back();
	perp_minBound = perp_bins.front(), perp_maxBound = perp_bins.back();
	para_cdf = gen_cdf(para_pdf,nbins,para_minBound,para_maxBound);  
	perp_cdf = gen_cdf(perp_pdf,nbins,perp_minBound,perp_maxBound);	
  }

  std::vector<double> operator()() {
    double para = cdf_sample(para_cdf,nbins,para_minBound,para_maxBound);
    double perp = cdf_sample(perp_cdf,nbins,perp_minBound,perp_maxBound);
    return ecliptic_to_equatorial(para,perp);
  }
	
 private:

  std::vector<double> perp_cdf, para_cdf;
  int nbins;
  double para_minBound, perp_minBound, para_maxBound, perp_maxBound;
  
  //function for sampling a random value from the inverse cdf (inversion is done implicitly)
  double cdf_sample(std::vector<double> cdf, const uint32_t &nbins, const double &minBound, const double maxBound)
  {
  	double r = drand48(); //figure out a way to include engine?
	std::vector<double>::iterator lwr = std::lower_bound(cdf.begin(), cdf.end(), r);
  	int off = std::max(0, (int)(lwr - cdf.begin()));
  	double t = (r -cdf[off])/(cdf[off + 1] - cdf[off]);
  	double x = (off + t)/(double)(nbins);
  	return minBound + (maxBound-minBound)*x;
  }
  
  //function for creating a cumulative distribution function from a probability distribution function
  std::vector<double> gen_cdf(std::vector<double> pdf, uint32_t &nbins, const double &minBound, const double maxBound)
  {
	std::vector<double> cdf;
	double dx = (maxBound - minBound)/nbins;
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
  std::vector<double> ecliptic_to_equatorial(float para, float perp)
  {
	double ecl_obl = 0.4090926; //radians
	double dec = asin(sin(ecl_obl)*sin(para/206265.0)*cos(perp/206265.0)+cos(ecl_obl)*sin(perp/206265.0));
	double ra = acos((cos(para/206265.0)*cos(perp/206265.0))/cos(dec));
	std::vector<double> v = {ra*206265.0,dec*2065265.0}; //now in arcsec/sec
	return v; 
  }
};	
	
} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_CUSTOM_DISTRIBUTION_HPP
