#include "dft.h"

#include <cmath>
#include <complex>
#include <vector>

namespace opendsp {

using namespace std::complex_literals;

std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>>& x)
{
    size_t N = x.size();
    std::vector<std::complex<double>> X(N);

    for (int k = 0; k < N; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            //X[k] += x[n] * (std::cos(2 * M_PI * n * k / N) - 1i * std::sin(2 * M_PI * n * k / N));
            X[k] += x[n] * std::exp(-1i * 2.0 * M_PI * static_cast<double>(n) * static_cast<double>(k) / static_cast<double>(N));
        }
    }

    return X;
}

std::vector<std::complex<double>> DFT(const std::vector<double>& x)
{
    size_t N = x.size();
    std::vector<std::complex<double>> X(N);

    // Take advantage of DFT symmetry when dealing with real input signals
    // Only the first N/2 + 1 outputs are unique
    for (int k = 0; k < N / 2 + 1; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            //X[k] += x[n] * (std::cos(2 * M_PI * n * k / N) - 1i * std::sin(2 * M_PI * n * k / N));
            X[k] += x[n] * std::exp(-1i * 2.0 * M_PI * static_cast<double>(n) * static_cast<double>(k) / static_cast<double>(N));
        }

        // X(N-k) = X(k)* for k = 1 -> N/2
        if (k != 0)
        {
            X[N - k] = std::conj(X[k]);
        }
    }

    return X;
}

} /* namespace opendsp */