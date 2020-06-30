#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <complex>

double now();
void writewav(const std::vector<double> &samples, const char *filename, int rate);
std::vector<double> readwav(const char *filename, int &rate_out);
void writetxt(std::vector<double> v, const char *filename);
std::vector<double> vrange(std::vector<double> v, int start, int len);
std::vector<double> raised_cosine(int n);
std::vector<double> convolve(const std::vector<double> &a, const std::vector<double> &b);
std::vector<std::complex<double>> cconvolve(const std::vector<std::complex<double>> &a, const std::vector<double> &b);
std::vector<double> thin(const std::vector<double> &a, int factor);
std::complex<double> goertzel(std::vector<double> v, int rate, int i0, int n, double hz);

typedef unsigned long ulong;
typedef unsigned int uint;

#endif
