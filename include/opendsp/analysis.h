#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

std::vector<double> Magnitude(std::vector<std::complex<double>> X);

std::vector<double> Phase(std::vector<std::complex<double>> X, double threshold, double amplitude);

} /* namespace opendsp */
