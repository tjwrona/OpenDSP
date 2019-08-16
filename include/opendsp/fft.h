#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

std::vector<std::complex<double>> FFT(const std::vector<double>& x);

} /* namespace opendsp */