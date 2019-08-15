#include "fft.h"

#include <cassert>
#include <complex>
#include <vector>

#include "dft.h"
#include "signal.h"

namespace opendsp
{

static Signal<std::complex<double>> FFTRadix2(const Signal<std::complex<double>>& x);
static bool IsPowerOf2(size_t x);
static size_t ReverseBits(const size_t x, const size_t bitCount);

//TODO: CLEAN THIS UP!!!
//      Also this requires some significant error handling
Signal<std::complex<double>> FFT(const Signal<double>& x)
{
    Signal<std::complex<double>> tmp(x.GetSampleRate(), x.GetLength() / 2);
    for (int i = 0; i < tmp.GetLength(); ++i)
    {
        tmp[i] = std::complex<double>(x[2 * i], x[2 * i + 1]);
    }

    Signal<std::complex<double>> Zn = FFTRadix2(tmp);

    Signal<std::complex<double>> Xn(Zn.GetSampleRate(), Zn.GetLength() + 1);
    Signal<std::complex<double>> Yn(Zn.GetSampleRate(), Zn.GetLength() + 1);

    Signal<std::complex<double>> X(x.GetSampleRate(), x.GetLength());

    size_t N = Zn.GetLength();
    Xn[0] = std::complex<double>(Zn[0].real(), 0);
    Yn[0] = std::complex<double>(Zn[0].imag(), 0);
    Xn[N] = std::complex<double>(Zn[N / 2].real(), 0);
    Yn[N] = std::complex<double>(Zn[N / 2].imag(), 0);

    for (int k = 1; k < N; ++k)
    {
        Xn[k] = std::complex<double>((Zn[k].real() + Zn[N - k].real()) / 2, (Zn[k].imag() - Zn[N - k].imag()) / 2);
        Yn[k] = std::complex<double>((Zn[N - k].imag() + Zn[k].imag()) / 2, (Zn[N - k].real() - Zn[k].real()) / 2);

        std::complex<double> W = std::polar(1.0, -2 * M_PI * k / x.GetLength());
        X[k] = Xn[k] + W * Yn[k];

        if (k != 0)
        {
            X[X.GetLength() - k] = std::conj(X[k]);
        }
    }

    return X;
}

Signal<std::complex<double>> FFT(const Signal<std::complex<double>>& x)
{
    if (IsPowerOf2(x.GetLength()))
    {
        return FFTRadix2(x);
    }
    else
    {
        return DFT(x);
    }
}

static Signal<std::complex<double>> FFTRadix2(const Signal<std::complex<double>>& x)
{
    assert(IsPowerOf2(x.GetLength()));

    size_t stages = static_cast<size_t>(log2(x.GetLength()));

    Signal<std::complex<double>> X(x.GetSampleRate(), x.GetLength());
    for (size_t i = 0; i < x.GetLength(); ++i)
    {
        X[i] = std::complex<double>(x[ReverseBits(i, stages)]);
    }

    for (size_t stage = 1; stage <= stages; ++stage)
    {
        size_t N = static_cast<size_t>(std::pow(2, stage));
        for (size_t k = 0; k < x.GetLength(); k += N)
        {
            for (size_t n = 0; n < N / 2; ++n)
            {
                //TODO: Avoid duplicate twiddle factor calcs
                //TODO: Also consider pre-calculating or cacheing them
                std::complex<double> W = std::polar(1.0, -2 * M_PI * (k + n) / N);
                std::complex<double> tmp = X[k + n];
                X[k + n] = tmp + W * X[k + n + N / 2];
                X[k + n + N / 2] = tmp - W * X[k + n + N / 2];
            }
        }
    }

    return X;
}

static bool IsPowerOf2(size_t x)
{
    return x && (!(x & (x - 1)));
}

static size_t ReverseBits(const size_t x, const size_t bitCount)
{
    size_t xReversed = 0;
    for (size_t i = 0; i < bitCount; ++i)
    {
        xReversed = (xReversed << 1U) | ((x >> i) & 1U);
    }

    return xReversed;
}

} /* namespace opendsp */