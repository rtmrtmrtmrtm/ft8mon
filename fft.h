#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

std::vector<std::complex<double>> one_fft(const std::vector<double> &samples,
                                          int i0, int block, const char *why);
std::vector<double> one_ifft(const std::vector<std::complex<double>> &bins, const char *why);

typedef std::vector< std::vector< std::complex<double> > > ffts_t;
ffts_t ffts(const std::vector<double> &samples, int i0, int block, const char *why);

std::vector<std::complex<double>> one_fft_c(const std::vector<double> &samples,
                                            int i0, int block, const char *why);
std::vector<std::complex<double>> one_ifft_cc(const std::vector<std::complex<double>> &bins, const char *why);

std::vector<std::complex<double>> analytic(const std::vector<double> &x, const char *why);

std::vector<double> hilbert_shift(const std::vector<double> &x,
                                  double hz0, double hz1, int rate);


#endif
