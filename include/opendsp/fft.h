#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

std::vector<std::complex<double>> FFT(const std::vector<double>& x);
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& x);

} /* namespace opendsp */