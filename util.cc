#include "sndfile.h"
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <complex>
#include "util.h"

double
now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void
writewav(const std::vector<double> &samples, const char *filename, int rate)
{
  double mx = 0;
  for(ulong i = 0; i < samples.size(); i++){
    if(fabs(samples[i]) > mx)
      mx = fabs(samples[i]);
  }
  std::vector<double> v(samples.size());
  for(ulong i = 0; i < samples.size(); i++){
    v[i] = (samples[i] / mx) * (2000.0 / 32767.0);
  }

  SF_INFO sf;
  sf.channels = 1;
  sf.samplerate = rate;
  sf.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  SNDFILE *f = sf_open(filename, SFM_WRITE, &sf);
  assert(f);
  sf_write_double(f, v.data(), v.size());
  sf_write_sync(f);
  sf_close(f);
}

std::vector<double>
readwav(const char *filename, int &rate_out)
{
  SF_INFO info;
  memset(&info, 0, sizeof(info));
  SNDFILE *sf = sf_open(filename, SFM_READ, &info);
  if(sf == 0){
    fprintf(stderr, "cannot open %s\n", filename);
    exit(1); // XXX
  }
  rate_out = info.samplerate;

  std::vector<double> out;

  while(1){
    double buf[512];
    int n = sf_read_double(sf, buf, 512);
    if(n <= 0)
      break;
    for(int i = 0; i < n; i++){
      out.push_back(buf[i]);
    }
  }

  sf_close(sf);

  return out;
}

void
writetxt(std::vector<double> v, const char *filename)
{
  FILE *fp = fopen(filename, "w");
  if(fp == 0){
    fprintf(stderr, "could not write %s\n", filename);
    exit(1);
  }
  for(ulong i = 0; i < v.size(); i++){
    fprintf(fp, "%f\n", v[i]);
  }
  fclose(fp);
}

//
// FIR coefficients for a "raised cosine" filter, which I believe
// is how BPSK31 does transmit filtering.
// https://ham.stackexchange.com/questions/7703/filtering-for-psk31-demodulation
//
// to be applied to the -1/+1 phases, not the
// actual signal.
//
// the coefficients sum to 1.0.
//
std::vector<double>
raised_cosine(int n)
{
  std::vector<double> v;
  double sum = 0.0;
  for(int i = 0; i < n; i++){
    double y = -cos(2 * M_PI * (i / (double) n));
    y += 1;
    y /= 2;
    v.push_back(y);
    sum += y;
  }
  for(int i = 0; i < n; i++){
    v[i] /= sum;
  }
  return v;
}

//
// a[] should be longer than b[].
//
std::vector<double>
convolve(const std::vector<double> &a, const std::vector<double> &b)
{
  int alen = a.size();
  int blen = b.size();
  std::vector<double> c;
  for(int i = 0; i < alen; i++){
    double sum = 0;
    for(int j = 0; j < blen; j++){
      int ii = i + blen/2 - j;
      if(ii >= 0 && ii < alen){
        sum += a[ii] * b[j];
      }
    }
    c.push_back(sum);
  }
  return c;
}


//
// a[] should be longer than b[].
//
std::vector<std::complex<double>>
cconvolve(const std::vector<std::complex<double>> &a, const std::vector<double> &b)
{
  int alen = a.size();
  int blen = b.size();
  std::vector<std::complex<double>> c;
  for(int i = 0; i < alen; i++){
    double rsum = 0;
    double isum = 0;
    for(int j = 0; j < blen; j++){
      int ii = i + blen/2 - j;
      if(ii >= 0 && ii < alen){
        rsum += a[ii].real() * b[j];
        isum += a[ii].imag() * b[j];
      }
    }
    c.push_back(std::complex(rsum, isum));
  }
  return c;
}

//
// every factor'th sample.
// does not low-pass filter.
//
std::vector<double>
thin(const std::vector<double> &a, int factor)
{
  assert((a.size() % factor) == 0);
  std::vector<double> out(a.size() / factor);
  for(int i = 0; i < out.size(); i++){
    out[i] = a[i * factor];
  }
  return out;
}

//
// https://www.embedded.com/the-goertzel-algorithm/
//
std::complex<double>
old_goertzel(std::vector<double> v, int rate, int i0, int n, double hz)
{
  // int k = 0.5 + (n * hz) / rate;
  double k = (n * hz) / rate;
  double w = (2.0 * M_PI / n) * k;
  double coeff = 2 * cos(w);

  double q1 = 0;
  double q2 = 0;
  for(int i = 0; i < n; i++){
    double q0 = coeff*q1 - q2 + v[i0+i];
    q2 = q1;
    q1 = q0;
  }

  double re = q1 - q2 * cos(w);
  double im = q2 * sin(w);

  return std::complex<double>(re, im);
}

//
// Goertzel Algorithm for a Non-integer Frequency Index, Rick Lyons
// https://www.dsprelated.com/showarticle/495.php
//
std::complex<double>
goertzel(std::vector<double> v, int rate, int i0, int n, double hz)
{
  //double radians_per_sample = (hz * 2 * M_PI) / rate;
  //double k = radians_per_sample * n;
  double bin_hz = rate / (double) n;
  double k = hz / bin_hz;

  double alpha = 2 * M_PI * k / n;
  double beta = 2 * M_PI * k * (n - 1.0) / n;

  double two_cos_alpha = 2 * cos(alpha);
  double a = cos(beta);
  double b = -sin(beta);
  double c = sin(alpha) * sin(beta) - cos(alpha)*cos(beta);
  double d = sin(2 * M_PI * k);

  double w1 = 0;
  double w2 = 0;

  for(int i = 0; i < n; i++){
    double w0 = v[i0+i] + two_cos_alpha * w1 - w2;
    w2 = w1;
    w1 = w0;
  }

  double re = w1*a + w2*c;
  double im = w1*b + w2*d;

  return std::complex<double>(re, im);
}
