#include "analysis.h"

#include <algorithm>
#include <complex>
#include <vector>

namespace opendsp {

std::vector<double> Magnitude(std::vector<std::complex<double>> X)
{
    size_t N = X.size();
    std::vector<double> X_magnitude(N);

    std::transform(X.begin(), X.end(), X_magnitude.begin(), [](const std::complex<double>& value) {
        return std::abs(value);
    });

    return X_magnitude;
}

std::vector<double> Phase(std::vector<std::complex<double>> X, double threshold, double amplitude)
{
    size_t N = X.size();
    std::vector<double> X_phase(N);

    std::transform(X.begin(), X.end(), X_phase.begin(), [threshold, amplitude](const std::complex<double>& value) {
        double level = 10.0 * std::log10(std::norm(value) / std::pow(amplitude, 2.0));
        return level > threshold ? std::arg(value) : 0.0;
    });

    return X_phase;
}

} /* namespace opendsp */