#include "fft.h"

#include <bitset>
#include <cassert>
#include <complex>
#include <vector>

namespace opendsp
{

static std::vector<std::complex<double>> FFTRadix2(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& W);
static bool IsPowerOf2(size_t x);

static constexpr size_t wordSize = sizeof(size_t) * 8;

std::vector<std::complex<double>> FFT(const std::vector<double>& x)
{
    const size_t N = x.size();

    // Radix2 FFT requires length of the input signal to be a power of 2
    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(N));
    const size_t NOver2 = N / 2;

    // Taking advantage of symmetry the FFT of a real signal can be computed
    // using a single N/2-point complex FFT. Split the input signal into its
    // even and odd components and load the data into a single complex vector.
    std::vector<std::complex<double>> x_p(NOver2);
    for (size_t n = 0; n < NOver2; ++n)
    {
        // x_p[n] = x[2n] + jx[2n + 1]
        const auto nTimes2 = n * 2;
        x_p[n] = std::complex<double>(x[nTimes2], x[nTimes2 + 1]);
    }

    // Pre-calculate twiddle factors
    std::vector<std::complex<double>> W(NOver2);
    std::vector<std::complex<double>> W_p(NOver2 / 2);
    const auto twiddleConstant = -2.0 * M_PI / N;
    for (size_t k = 0; k < NOver2; ++k)
    {
        W[k] = std::polar(1.0, k * twiddleConstant);

        // The N/2-point complex DFT uses only the even twiddle factors
        if (k % 2 == 0)
        {
            W_p[k / 2] = W[k];
        }
    }

    // Perform the N/2-point complex FFT
    std::vector<std::complex<double>> X_p = FFTRadix2(x_p, W_p);

    // Extract the N-point FFT of the real signal from the results 
    std::vector<std::complex<double>> X(N);
    X[0] = X_p[0].real() + X_p[0].imag();
    for (size_t k = 1; k < NOver2; ++k)
    {
        const auto l = NOver2 - k;

        // Extract the FFT of the even components
        const auto A = std::complex<double>(
            (X_p[k].real() + X_p[l].real()) / 2,
            (X_p[k].imag() - X_p[l].imag()) / 2);

        // Extract the FFT of the odd components
        const auto B = std::complex<double>(
            (X_p[l].imag() + X_p[k].imag()) / 2,
            (X_p[l].real() - X_p[k].real()) / 2);

        // Sum the results and take advantage of symmetry
        const auto tmp = W[k] * B;
        X[k] = A + tmp;
        X[k + NOver2] = A - tmp;
    }

    return X;
}

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& x)
{
    const size_t N = x.size();

    // Radix2 FFT requires length of the input signal to be a power of 2
    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(N));
    const size_t NOver2 = N / 2;

    // Pre-calculate twiddle factors
    std::vector<std::complex<double>> W(NOver2);
    const auto twiddleConstant = -2.0 * M_PI / N;
    for (size_t k = 0; k < NOver2; ++k)
    {
        W[k] = std::polar(1.0, k * twiddleConstant);
    }

    return FFTRadix2(x, W);
}

static std::vector<std::complex<double>> FFTRadix2(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& W)
{
    const size_t N = x.size();

    // Radix2 FFT requires length of the input signal to be a power of 2
    assert(IsPowerOf2(N));
    const size_t NOver2 = N / 2;

    // Calculate how many stages the FFT must compute
    const size_t stages = static_cast<size_t>(log2(N));

    // Pre-load the output vector with the input data using bit-reversed indexes
    std::vector<std::complex<double>> X(N);
    size_t nReversed = 0;
    for (size_t n = 0; n < N; ++n)
    {
        X[n] = x[nReversed];
        const size_t change = n ^ (n + 1);
        const std::bitset<wordSize> bits(~change);
        nReversed ^= change << (bits.count() - (wordSize - stages));
    }

    // Calculate the FFT one stage at a time and sum the results
    size_t NPreviousStage = 1;
    size_t NStage = 2;
    size_t WOffset = NOver2;
    for (size_t stage = 1; stage <= stages; ++stage)
    {
        for (size_t k = 0; k < N; k += NStage)
        {
            for (size_t n = 0; n < NPreviousStage; ++n)
            {
                const auto index1 = k + n;
                const auto index2 = index1 + NPreviousStage;
                const auto tmp1 = X[index1];
                const auto tmp2 = W[n * WOffset] * X[index2];
                X[index1] = tmp1 + tmp2;
                X[index2] = tmp1 - tmp2;
            }
        }
        NPreviousStage = NStage;
        NStage *= 2;
        WOffset /= 2;
    }

    return X;
}

// Returns true if x is a power of 2
static bool IsPowerOf2(size_t x)
{
    return x && (!(x & (x - 1)));
}

} /* namespace opendsp */