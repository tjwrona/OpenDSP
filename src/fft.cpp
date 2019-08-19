#include "fft.h"

#include <cassert>
#include <complex>
#include <vector>

namespace opendsp
{

using std::complex;
using std::polar;
using std::vector;

static vector<complex<double>> FFTStockham(const vector<complex<double>>& x, const vector<complex<double>>& W);
static bool IsPowerOf2(size_t x);

vector<complex<double>> FFT(const vector<double>& x)
{
    const size_t N = x.size();

    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(N));

    const size_t NOver2 = N / 2;

    // Taking advantage of symmetry the FFT of a real signal can be computed
    // using a single N/2-point complex FFT. Split the input signal into its
    // even and odd components and load the data into a single complex vector.
    vector<complex<double>> x_p(NOver2);
    for (size_t n = 0; n < NOver2; ++n)
    {
        // x_p[n] = x[2n] + jx[2n + 1]
        const auto nTimes2 = n * 2;
        x_p[n] = complex<double>(x[nTimes2], x[nTimes2 + 1]);
    }

    // Pre-calculate twiddle factors
    vector<complex<double>> W(NOver2);
    vector<complex<double>> W_p(NOver2 / 2);
    const auto twiddleConstant = -2.0 * M_PI / N;
    for (size_t k = 0; k < NOver2; ++k)
    {
        W[k] = polar(1.0, k * twiddleConstant);

        // The N/2-point complex FFT uses only the even twiddle factors
        if (k % 2 == 0)
        {
            W_p[k / 2] = W[k];
        }
    }

    // Perform the N/2-point complex FFT
    vector<complex<double>> X_p = FFTStockham(x_p, W_p);

    // Extract the N-point FFT of the real signal from the results 
    vector<complex<double>> X(N);
    X[0] = X_p[0].real() + X_p[0].imag();
    for (size_t k = 1; k < NOver2; ++k)
    {
        const auto NOver2MinusK = NOver2 - k;

        // Extract the FFT of the even components
        const auto A = complex<double>(
            (X_p[k].real() + X_p[NOver2MinusK].real()) / 2,
            (X_p[k].imag() - X_p[NOver2MinusK].imag()) / 2);

        // Extract the FFT of the odd components
        const auto B = complex<double>(
            (X_p[NOver2MinusK].imag() + X_p[k].imag()) / 2,
            (X_p[NOver2MinusK].real() - X_p[k].real()) / 2);

        // Sum the results
        X[k] = A + W[k] * B;
        X[N - k] = conj(X[k]); // (>*.*)> symmetry! <(*.*<)
    }

    return X;
}

vector<complex<double>> FFT(const vector<complex<double>>& x)
{
    const size_t N = x.size();

    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(N));

    const size_t NOver2 = N / 2;

    // Pre-calculate twiddle factors
    vector<complex<double>> W(NOver2);
    const auto twiddleConstant = -2.0 * M_PI / N;
    for (size_t k = 0; k < NOver2; ++k)
    {
        W[k] = polar(1.0, k * twiddleConstant);
    }

    return FFTStockham(x, W);
}

static vector<complex<double>> FFTStockham(const vector<complex<double>>& x, const vector<complex<double>>& W)
{
    const size_t N = x.size();

    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(N));

    const size_t NOver2 = N / 2;

    vector<complex<double>> a(x);
    vector<complex<double>> b(N);

    size_t WStride = NOver2;
    for (size_t stride = 1; stride < N; stride *= 2)
    {
        for (size_t k = 0; k < NOver2; k += stride)
        {
            const auto kTimes2 = k * 2;
            for (size_t n = 0; n < stride; ++n)
            {
                const auto aIndex1 = n + k;
                const auto aIndex2 = aIndex1 + NOver2;

                const auto bIndex1 = n + kTimes2;
                const auto bIndex2 = bIndex1 + stride;

                const auto tmp1 = a[aIndex1];
                const auto tmp2 = W[n * WStride] * a[aIndex2];

                b[bIndex1] = tmp1 + tmp2;
                b[bIndex2] = tmp1 - tmp2; // (>*.*)> symmetry! <(*.*<)
            }
        }

        WStride /= 2;
        a.swap(b);
    }

    return a;
}

// Returns true if x is a power of 2
static bool IsPowerOf2(size_t x)
{
    return x && (!(x & (x - 1)));
}

} /* namespace opendsp */