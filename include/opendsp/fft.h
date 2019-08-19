#pragma once

#include <complex>
#include <vector>

namespace opendsp
{

/**
 * @brief Calculates the fast Fourier transform of a real input signal
 *
 * Currently only a Radix-2 FFT has been implemented. The Radix-2 algorithm
 * is one of the fastest, but it has the restriction that the input vector must
 * have a length that is a power of 2.
 *
 * @param[in] x - Vector of real input data (whose length is a power of 2)
 * @return Vector of complex values representing the FFT of the input data
 */
std::vector<std::complex<double>> FFT(const std::vector<double>& x);

/**
 * @brief Calculates the fast Fourier transform of a complex inpout signal
 *
 * Currently only a Radix-2 FFT has been implemented. The Radix-2 algorithm
 * is one of the fastest, but it has the restriction that the input vector must
 * have a length that is a power of 2.
 *
 * @param[in] x - Vector of complex input data (whose length is a power of 2)
 * @return Vector of complex values representing the FFT of the input data
 */
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& x);

} /* namespace opendsp */