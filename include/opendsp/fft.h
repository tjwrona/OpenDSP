#pragma once

#include <complex>

#include "signal.h"

namespace opendsp
{

Signal<std::complex<double>> FFT(const Signal<double>& x);

} /* namespace opendsp */