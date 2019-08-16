#include <cmath>
#include <gtest/gtest.h>
#include <opendsp/dft.h>
#include <opendsp/fft.h>

//TODO: Remove
#include <chrono>

using namespace opendsp;

//TODO: Implement real tests
TEST(FFTTest, Test1)
{
    int N = 65536;
    std::vector<double> x(N);

    int f_s = 8000;
    double t_s = 1.0 / f_s;

    for (int n = 0; n < N; ++n)
    {
        x[n] = std::sin(2 * M_PI * 1000 * n * t_s) + 0.5 * std::sin(2 * M_PI * 2000 * n * t_s + 3 * M_PI / 4);
    }

    std::cout << "start" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto X = FFT(x);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "end - duration: " << duration.count() << " microseconds." << std::endl;

    //auto X2 = DFT(x);
    //auto X_magnitude = Magnitude(X);
    //auto X_power = Power(X);
    //auto X_phase = Phase(X, -120, 1.0);
}