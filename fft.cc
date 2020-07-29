#include "fft.h"
#include <fftw3.h>
#include <mutex>
#include <unistd.h>
#include <assert.h>
#include "util.h"

int fftw_type = FFTW_MEASURE;

// a cached fftw plan, for both of:
// fftw_plan_dft_r2c_1d(n, m_in, m_out, FFTW_ESTIMATE);
// fftw_plan_dft_c2r_1d(n, m_in, m_out, FFTW_ESTIMATE);
class Plan {
public:
  int n_;

  //
  // real -> complex
  //
  fftw_complex *c_; // (n_ / 2) + 1 of these
  double *r_; // n_ of these
  fftw_plan fwd_; // forward plan
  fftw_plan rev_; // reverse plan

  //
  // complex -> complex
  //
  fftw_complex *cc1_; // n
  fftw_complex *cc2_; // n
  fftw_plan cfwd_; // forward plan
  fftw_plan crev_; // reverse plan
};

static std::mutex plansmu;
static std::vector<Plan> plans;
static int plan_master_pid;

void
get_plan(int n, Plan &p)
{
  // cache fftw plans in the parent process,
  // so they will already be there for fork()ed children.
  
  plansmu.lock();
  
  if(plan_master_pid < 0){
    plan_master_pid = getpid();
  }
    
  for(ulong i = 0; i < plans.size(); i++){
    if(plans[i].n_ == n){
      p = plans[i];
      plansmu.unlock();
      return;
    }
  }
  
  //
  // real -> complex
  //
  
  p.n_ = n;
  p.r_ = (double*) fftw_malloc(n * sizeof(double));
  assert(p.r_);
  p.c_ = (fftw_complex*) fftw_malloc(((n/2)+1) * sizeof(fftw_complex));
  assert(p.c_);
  
  // FFTW_ESTIMATE
  // FFTW_MEASURE
  // FFTW_PATIENT
  // FFTW_EXHAUSTIVE
  int type = fftw_type;
  if(getpid() != plan_master_pid){
    type = FFTW_ESTIMATE;
  }
  p.fwd_ = fftw_plan_dft_r2c_1d(n, p.r_, p.c_, type);
  assert(p.fwd_);
  p.rev_ = fftw_plan_dft_c2r_1d(n, p.c_, p.r_, type);
  assert(p.rev_);
  
  //
  // complex -> complex
  //
  p.cc1_ = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
  assert(p.cc1_);
  p.cc2_ = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
  assert(p.cc2_);
  p.cfwd_ = fftw_plan_dft_1d(n, p.cc1_, p.cc2_, FFTW_FORWARD, type);
  assert(p.cfwd_);
  p.crev_ = fftw_plan_dft_1d(n, p.cc2_, p.cc1_, FFTW_BACKWARD, type);
  assert(p.crev_);

  plans.push_back(p);

  plansmu.unlock();
}

//
// do just one FFT on samples[i0..i0+block]
// real inputs, complex outputs.
// output has (block / 2) + 1 points.
//
std::vector<std::complex<double>>
one_fft(const std::vector<double> &samples, int i0, int block)
{
  assert(i0 >= 0);
  assert(block > 1);
  
  int nsamples = samples.size();
  int nbins = (block / 2) + 1;

  Plan p;
  get_plan(block, p);
  fftw_plan m_plan = p.fwd_;

  double *m_in = (double *) fftw_malloc(sizeof(double) * p.n_);
  assert(m_in);
  fftw_complex *m_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *
                                                     ((p.n_ / 2) + 1));
  assert(m_out);

  for(int i = 0; i < block; i++){
    if(i0 + i < nsamples){
      m_in[i] = samples[i0 + i];
    } else {
      m_in[i] = 0;
    }
  }

  fftw_execute_dft_r2c(m_plan, m_in, m_out);

  std::vector<std::complex<double>> out(nbins);

  for(int bi = 0; bi < nbins; bi++){
    double re = m_out[bi][0];
    double im = m_out[bi][1];
    std::complex<double> c(re, im);
    out[bi] = c;
  }

  fftw_free(m_in);
  fftw_free(m_out);

  return out;
}

//
// do a full set of FFTs, one per symbol-time.
// bins[time][frequency]
//
ffts_t
ffts(const std::vector<double> &samples, int i0, int block)
{
  assert(i0 >= 0);
  assert(block > 1 && (block % 2) == 0);
  
  int nsamples = samples.size();
  int nbins = (block / 2) + 1;
  int nblocks = (nsamples - i0) / block;
  ffts_t bins(nblocks);
  for(int si = 0; si < nblocks; si++){
    bins[si].resize(nbins);
  }

  Plan p;
  get_plan(block, p);
  fftw_plan m_plan = p.fwd_;

  // allocate our own b/c using p.m_in and p.m_out isn't thread-safe.
  double *m_in = (double *) fftw_malloc(sizeof(double) * p.n_);
  fftw_complex *m_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *
                                                     ((p.n_ / 2) + 1));
  assert(m_in && m_out);

  // double *m_in = p.r_;
  // fftw_complex *m_out = p.c_;

  for(int si = 0; si < nblocks; si++){
    int off = i0 + si * block;
    for(int i = 0; i < block; i++){
      if(off + i < nsamples){
        double x = samples[off + i];
        m_in[i] = x;
      } else {
        m_in[i] = 0;
      }
    }

    fftw_execute_dft_r2c(m_plan, m_in, m_out);

    for(int bi = 0; bi < nbins; bi++){
      double re = m_out[bi][0];
      double im = m_out[bi][1];
      std::complex<double> c(re, im);
      bins[si][bi] = c;
    }
  }

  fftw_free(m_in);
  fftw_free(m_out);

  return bins;
}

//
// do just one FFT on samples[i0..i0+block]
// real inputs, complex outputs.
// output has block points.
//
std::vector<std::complex<double>>
one_fft_c(const std::vector<double> &samples, int i0, int block)
{
  assert(i0 >= 0);
  assert(block > 1);
  
  int nsamples = samples.size();

  Plan p;
  get_plan(block, p);
  fftw_plan m_plan = p.cfwd_;

  fftw_complex *m_in  = (fftw_complex*) fftw_malloc(block * sizeof(fftw_complex));
  fftw_complex *m_out = (fftw_complex*) fftw_malloc(block * sizeof(fftw_complex));
  assert(m_in && m_out);

  for(int i = 0; i < block; i++){
    if(i0 + i < nsamples){
      m_in[i][0] = samples[i0 + i]; // real
    } else {
      m_in[i][0] = 0;
    }
    m_in[i][1] = 0; // imaginary
  }

  fftw_execute_dft(m_plan, m_in, m_out);

  std::vector<std::complex<double>> out(block);

  double norm = 1.0 / sqrt(block);
  for(int bi = 0; bi < block; bi++){
    double re = m_out[bi][0];
    double im = m_out[bi][1];
    std::complex<double> c(re, im);
    c *= norm;
    out[bi] = c;
  }
    
  fftw_free(m_in);
  fftw_free(m_out);

  return out;
}

std::vector<std::complex<double>>
one_ifft_cc(const std::vector<std::complex<double>> &bins)
{
  int block = bins.size();

  Plan p;
  get_plan(block, p);
  fftw_plan m_plan = p.crev_;

  fftw_complex *m_in = (fftw_complex*) fftw_malloc(block * sizeof(fftw_complex));
  fftw_complex *m_out = (fftw_complex *) fftw_malloc(block * sizeof(fftw_complex));
  assert(m_in && m_out);

  for(int bi = 0; bi < block; bi++){
    double re = bins[bi].real();
    double im = bins[bi].imag();
    m_in[bi][0] = re;
    m_in[bi][1] = im;
  }

  fftw_execute_dft(m_plan, m_in, m_out);

  std::vector<std::complex<double>> out(block);
  double norm = 1.0 / sqrt(block);
  for(int i = 0; i < block; i++){
    double re = m_out[i][0];
    double im = m_out[i][1];
    std::complex<double> c(re, im);
    c *= norm;
    out[i] = c;
  }

  fftw_free(m_in);
  fftw_free(m_out);

  return out;
}

std::vector<double>
one_ifft(const std::vector<std::complex<double>> &bins)
{
  int nbins = bins.size();
  int block = (nbins - 1) * 2;

  Plan p;
  get_plan(block, p);
  fftw_plan m_plan = p.rev_;

  fftw_complex *m_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *
                                                    ((p.n_ / 2) + 1));
  double *m_out = (double *) fftw_malloc(sizeof(double) * p.n_);

  for(int bi = 0; bi < nbins; bi++){
    double re = bins[bi].real();
    double im = bins[bi].imag();
    m_in[bi][0] = re;
    m_in[bi][1] = im;
  }

  fftw_execute_dft_c2r(m_plan, m_in, m_out);

  std::vector<double> out(block);
  for(int i = 0; i < block; i++){
    out[i] = m_out[i];
  }

  fftw_free(m_in);
  fftw_free(m_out);

  return out;
}

//
// return the analytic signal for signal x,
// just like scipy.signal.hilbert(), from which
// this code is copied.
//
// the return value is x + iy, where y is the hilbert transform of x.
//
std::vector<std::complex<double>>
analytic(const std::vector<double> &x)
{
  ulong n = x.size();

  std::vector<std::complex<double>> y = one_fft_c(x, 0, n);
  assert(y.size() == n);

  // leave y[0] alone.
  // double the first (positive) half of the spectrum.
  // zero out the second (negative) half of the spectrum.
  // y[n/2] is the nyquist bucket if n is even; leave it alone.
  if((n % 2) == 0){
    for(ulong i = 1; i < n/2; i++)
      y[i] *= 2;
    for(ulong i = n/2+1; i < n; i++)
      y[i] = 0;
  } else {
    for(ulong i = 1; i < (n+1)/2; i++)
      y[i] *= 2;
    for(ulong i = (n+1)/2; i < n; i++)
      y[i] = 0;
  }
      
  std::vector<std::complex<double>> z = one_ifft_cc(y);

  return z;
}

//
// general-purpose shift x in frequency by hz.
// uses hilbert transform to avoid sidebands.
// but it does wrap around at 0 hz and the nyquist frequency.
//
// note analytic() does an FFT over the whole signal, which
// is expensive, and often re-used, but it turns out it
// isn't a big factor in overall run-time.
//
// like weakutil.py's freq_shift().
//
std::vector<double>
hilbert_shift(const std::vector<double> &x, double hz0, double hz1, int rate)
{
  // y = scipy.signal.hilbert(x)
  std::vector<std::complex<double>> y = analytic(x);
  assert(y.size() == x.size());

  double dt = 1.0 / rate;
  int n = x.size();

  std::vector<double> ret(n);
  
  for(int i = 0; i < n; i++){
    // complex "local oscillator" at hz.
    double hz = hz0 + (i / (double)n) * (hz1 - hz0);
    std::complex<double> lo = std::exp(std::complex<double>(0.0, 2 * M_PI * hz * dt * i));
    ret[i] = (lo * y[i]).real();
  }

  return ret;
}
