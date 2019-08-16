#include <cmath>
#include <gtest/gtest.h>
#include <opendsp/dft.h>
#include <opendsp/fft.h>

using namespace opendsp;

TEST(FFT, FFT)
{
    size_t N = 8;
    std::vector<double> x(N);

    int f_s = 8000;
    double t_s = 1.0 / f_s;

    for (size_t n = 0; n < N; ++n)
    {
        x[n] = std::sin(2 * M_PI * 1000 * n * t_s)
             + 0.5 * std::sin(2 * M_PI * 2000 * n * t_s + 3 * M_PI / 4);
    }

    auto X = FFT(x);
    auto X2 = DFT(x);

    ASSERT_EQ(X.size(), N);
    ASSERT_EQ(X2.size(), N);
    
    double epsilon = std::pow(10.0, -std::numeric_limits<double>::digits10);

    for (size_t k = 0; k < N; ++k)
    {
        EXPECT_NEAR(X[k].real(), X2[k].real(), epsilon);
        EXPECT_NEAR(X[k].imag(), X2[k].imag(), epsilon);
    }
}