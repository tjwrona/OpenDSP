#include "dft.h"

#include <cmath>
#include <complex>
#include <vector>

#include "signal.h"

namespace opendsp
{

using namespace std::complex_literals;

Signal<std::complex<double>> DFT(const Signal<double>& x)
{
    Signal<std::complex<double>> X(x.GetSampleRate(), x.GetLength());

    // Take advantage of DFT symmetry when dealing with real input signals
    // Only the first N/2 + 1 outputs are unique
    size_t N = x.GetLength();
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

Signal<std::complex<double>> DFT(const Signal<std::complex<double>>& x)
{
    Signal<std::complex<double>> X(x.GetSampleRate(), x.GetLength());

    size_t N = x.GetLength();
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

} /* namespace opendsp */