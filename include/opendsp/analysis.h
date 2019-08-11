#ifndef OPENDSP_ANALYSIS_H
#define OPENDSP_ANALYSIS_H

#include <complex>
#include <vector>

namespace opendsp {

std::vector<double> Magnitude(std::vector<std::complex<double>> X);

std::vector<double> Phase(std::vector<std::complex<double>> X, double SNR, double amplitude);

} /* namespace opendsp */

#endif