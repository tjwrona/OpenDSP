#include <cmath>
#include <gtest/gtest.h>
#include <opendsp/analysis.h>
#include <opendsp/dft.h>

using namespace opendsp;

TEST(Test1, Test1)
{
    int N = 8;
    std::vector<double> x(N);

    int f_s = 8000;
    double t_s = 1.0 / f_s;

    for (int n = 0; n < N; ++n)
    {
        x[n] = std::sin(2 * M_PI * 1000 * n * t_s) + 0.5 * std::sin(2 * M_PI * 2000 * n * t_s + 3 * M_PI / 4);
    }

    auto X = DFT(x);
    auto X_magnitude = Magnitude(X);
    //auto X_power = Power(X);
    auto X_phase = Phase(X, -120, 1.0);
}