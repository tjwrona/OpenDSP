#include "analysis.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace opendsp
{

double Magnitude(const double value)
{
    return std::abs(value);
}

double Magnitude(const std::complex<double>& value)
{
    return std::abs(value);
}

std::vector<double> Magnitude(const std::vector<double>& x)
{
    std::vector<double> x_magnitude(x.size());
    std::transform(x.begin(), x.end(), x_magnitude.begin(), [](const double value) {
        return Magnitude(value);
        });

    return x_magnitude;
}

std::vector<double> Magnitude(const std::vector<std::complex<double>>& x)
{
    std::vector<double> x_magnitude(x.size());
    std::transform(x.begin(), x.end(), x_magnitude.begin(), [](const std::complex<double>& value) {
        return Magnitude(value);
        });

    return x_magnitude;
}

double Phase(const std::complex<double>& value)
{
    return std::arg(value);
}

std::vector<double> Phase(const std::vector<std::complex<double>>& x)
{
    std::vector<double> x_phase(x.size());
    std::transform(x.begin(), x.end(), x_phase.begin(), [](const std::complex<double>& value) {
        return Phase(value);
        });

    return x_phase;
}

double Power(const double value)
{
    return std::norm(value);
}

double Power(const std::complex<double>& value)
{
    return std::norm(value);
}

std::vector<double> Power(const std::vector<double>& x)
{
    std::vector<double> x_power(x.size());
    std::transform(x.begin(), x.end(), x_power.begin(), [](const double value) {
        return Power(value);
        });

    return x_power;
}

std::vector<double> Power(const std::vector<std::complex<double>>& x)
{
    std::vector<double> x_power(x.size());
    std::transform(x.begin(), x.end(), x_power.begin(), [](const std::complex<double>& value) {
        return Power(value);
        });

    return x_power;
}

} /* namespace opendsp */