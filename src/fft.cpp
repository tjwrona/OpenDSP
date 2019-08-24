#include "fft.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <vector>

namespace opendsp
{

using std::complex;
using std::conj;
using std::polar;
using std::transform;
using std::vector;

static bool IsPowerOf2(const size_t value);
static vector<complex<double>> StockhamFFT(const vector<double>& x);
static vector<complex<double>> StockhamFFT(const vector<complex<double>>& x);
static vector<complex<double>> StockhamFFT(const vector<complex<double>>& x, const vector<complex<double>>& W);

vector<complex<double>> FFT(const vector<double>& x)
{
    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(x.size()));
    return StockhamFFT(x);
}

vector<complex<double>> FFT(const vector<complex<double>>& x)
{
    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(x.size()));
    return StockhamFFT(x);
}

vector<complex<double>> IFFT(const vector<complex<double>>& X)
{
    // TODO: Implement other algorithms for when N is not a power of 2
    assert(IsPowerOf2(X.size()));

    vector<complex<double>> X_p(X.size());
    transform(X.begin(), X.end(), X_p.begin(), [](const complex<double>& value) {
        return conj(value);
        });

    const auto x_p = StockhamFFT(X_p);

    vector<complex<double>> x(x_p.size());
    transform(x_p.begin(), x_p.end(), x.begin(), [](const complex<double>& value) {
        return conj(value);
        });
}

// Returns true if x is a power of 2
static bool IsPowerOf2(const size_t value)
{
    return value && (!(value & (value - 1)));
}

vector<complex<double>> StockhamFFT(const vector<double>& x)
{
    assert(IsPowerOf2(x.size()));
    const size_t N = x.size();
    const size_t NOver2 = N / 2;

    // Taking advantage of symmetry the FFT of a real signal can be computed
    // using a single N/2-point complex FFT. Split the input signal into its
    // even and odd components and load the data into a single complex vector.
    vector<complex<double>> x_p(NOver2);
    for (size_t n = 0; n < NOver2; ++n)
    {
        const auto nTimes2 = n * 2;

        // x_p[n] = x[2n] + jx[2n + 1]
        x_p[n] = complex<double>(x[nTimes2], x[nTimes2 + 1]);
    }

    // Pre-calculate twiddle factors
    vector<complex<double>> W(NOver2);
    vector<complex<double>> W_p(NOver2 / 2);
    const auto omega = 2.0 * M_PI / N;
    for (size_t n = 0; n < NOver2; ++n)
    {
        W[n] = polar(1.0, -omega * n);

        // The N/2-point complex FFT uses only the even twiddle factors
        if (n % 2 == 0)
        {
            W_p[n / 2] = W[n];
        }
    }

    // Perform the N/2-point complex FFT
    vector<complex<double>> X_p = StockhamFFT(x_p, W_p);

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

static vector<complex<double>> StockhamFFT(const vector<complex<double>>& x)
{
    assert(IsPowerOf2(x.size()));
    const size_t N = x.size();
    const size_t NOver2 = N / 2;

    // Pre-calculate twiddle factors
    vector<complex<double>> W(NOver2);
    const auto omega = 2.0 * M_PI / N;
    for (size_t n = 0; n < NOver2; ++n)
    {
        W[n] = polar(1.0, -omega * n);
    }

    return StockhamFFT(x, W);
}

static vector<complex<double>> StockhamFFT(const vector<complex<double>>& x, const vector<complex<double>>& W)
{
    assert(IsPowerOf2(x.size()));
    const size_t N = x.size();
    const size_t NOver2 = N / 2;

    // The Stockham algorithm requires one vector for input/output data and
    // another as a temporary workspace
    vector<complex<double>> a(x);
    vector<complex<double>> b(N);

    // Set the spacing between twiddle factors used at the first stage
    size_t WStride = NOver2;

    // Loop through each stage of the FFT
    for (size_t stride = 1; stride < N; stride *= 2)
    {
        // Loop through the individual FFTs of each stage
        for (size_t k = 0; k < NOver2; k += stride)
        {
            const auto kTimes2 = k * 2;

            // Perform each individual FFT
            for (size_t n = 0; n < stride; ++n)
            {
                // Calculate the input indexes
                const auto aIndex1 = n + k;
                const auto aIndex2 = aIndex1 + NOver2;

                // Calculate the output indexes
                const auto bIndex1 = n + kTimes2;
                const auto bIndex2 = bIndex1 + stride;

                // Perform the FFT
                const auto tmp1 = a[aIndex1];
                const auto tmp2 = W[n * WStride] * a[aIndex2];

                // Sum the results
                b[bIndex1] = tmp1 + tmp2;
                b[bIndex2] = tmp1 - tmp2; // (>*.*)> symmetry! <(*.*<)
            }
        }

        // Spacing between twiddle factors is half for the next stage
        WStride /= 2;

        // Swap the data (output of this stage is input of the next)
        a.swap(b);
    }

    return a;
}

} /* namespace opendsp */