#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

double Magnitude(const double value);
double Magnitude(const std::complex<double>& value);
std::vector<double> Magnitude(const std::vector<double>& x);
std::vector<double> Magnitude(const std::vector<std::complex<double>>& x);

// Phase calculation doesn't make sense for real values
double Phase(const std::complex<double>& value);
std::vector<double> Phase(const std::vector<std::complex<double>>& x);

double Power(const double value);
double Power(const std::complex<double>& value);
std::vector<double> Power(const std::vector<double>& x);
std::vector<double> Power(const std::vector<std::complex<double>>& x);

} /* namespace opendsp */
