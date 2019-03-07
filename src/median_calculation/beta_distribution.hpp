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

#ifndef UMAP_APPS_MEDIAN_CALCULATION_BETA_DISTRIBUTION_HPP
#define UMAP_APPS_MEDIAN_CALCULATION_BETA_DISTRIBUTION_HPP

#include <random>

namespace median {

class beta_distribution {
 public:
  beta_distribution(double a, double b)
      : m_x_gamma(a, 1.0),
        m_y_gamma(b, 1.0) {}

  // () operator written to return two values for syntax help in the run_random_vector code (not often used)
  template <typename rnd_engine>
  std::vector<double> operator()(rnd_engine &engine) {
    double x1 = m_x_gamma(engine);
    double y1 = m_y_gamma(engine);
    double x2 = m_x_gamma(engine);
    double y2 = m_x_gamma(engine);
    std::vector<double> result = {(x1 / (x1 + y1)), (x2 / (x2 + y2))};
    return result;
  }

 private:
  std::gamma_distribution<> m_x_gamma;
  std::gamma_distribution<> m_y_gamma;
};

} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_BETA_DISTRIBUTION_HPP
