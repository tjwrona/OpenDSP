#include "analysis.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace opendsp
{

std::vector<double> Magnitude(std::vector<std::complex<double>> X)
{
    size_t N = X.size();
    std::vector<double> X_magnitude(N);

    std::transform(X.begin(), X.end(), X_magnitude.begin(), [](const std::complex<double>& value) {
        return std::abs(value);
        });

    return X_magnitude;
}

std::vector<double> Phase(std::vector<std::complex<double>> X)
{
    size_t N = X.size();
    std::vector<double> X_phase(N);

    std::transform(X.begin(), X.end(), X_phase.begin(), [](const std::complex<double>& value) {
        return std::arg(value);
        });

    return X_phase;
}

std::vector<double> Power(std::vector<std::complex<double>> X)
{
    size_t N = X.size();
    std::vector<double> X_power(N);

    std::transform(X.begin(), X.end(), X_power.begin(), [](const std::complex<double>& value) {
        return std::norm(value);
        });

    return X_power;
}

} /* namespace opendsp */