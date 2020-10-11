
//
// The core of an FT8 decoder in C++.
//
// Many ideas and protocol details borrowed from Franke
// and Taylor's WSJT-X code.
//
// Robert Morris, AB1HL
//

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <complex>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <random>
#include <functional>
#include <map>
#include "util.h"
#include "fft.h"
#include "ft8.h"

// 1920-point FFT at 12000 samples/second
// 6.25 Hz spacing, 0.16 seconds/symbol
// encode chain:
//   77 bits
//   append 14 bits CRC (for 91 bits)
//   LDPC(174,91) yields 174 bits
//   that's 58 3-bit FSK-8 symbols
//   gray code each 3 bits
//   insert three 7-symbol Costas sync arrays
//     at symbol #s 0, 36, 72 of final signal
//   thus: 79 FSK-8 symbols
// total transmission time is 12.64 seconds

// tunable parameters
int nthreads = 2; // number of parallel threads, for multi-core
int npasses = 3;  // number of spectral subtraction passes
int ldpc_iters = 25; // how hard LDPC decoding should work
int snr_win = 7; // averaging window, in symbols, for SNR conversion
int snr_how = 0; // technique to measure "N" for SNR. 0 means median of the 8 tones.
int soft_ranges = 1;
int best_in_noise = 1;
double shoulder200 = 10; // for 200 sps bandpass filter
double shoulder200_extra = 0.0; // for bandpass filter
double second_hz_win = 3.5; // +/- hz
int second_hz_fracs = 10; // divide total window into this many pieces
double second_off_win = 0.5; // +/- search window in symbol-times
int second_off_fracs = 10;
double third_hz_inc = 6.25 / 32;
double third_hz_win = 1;
int third_off_inc = 1;
int third_off_win = 16;
double log_tail = 0.1;
double log_rate = 8.0;
int problt_how = 3;
int use_apriori = 1;
int use_hints = 2; // 1 means use all hints, 2 means just CQ hints
int win_type = 1;
int osd_depth = 6; // don't increase beyond 6, produces too much garbage
int osd_ldpc_thresh = 70; // demand this many correct LDPC parity bits before OSD
int ncoarse = 1; // number of offsets per hz produced by coarse()
int ncoarse_blocks = 1;
double tminus = 2.2; // start looking at 0.5 - tminus seconds
double tplus = 2.4;
int coarse_off_fracs = 4;
int coarse_hz_fracs = 8;
double already_hz = 27;
double overlap = 20;
int overlap_edges = 0;
int sub_amp_win = 1;
int sub_phase_win = 0;
double nyquist = 0.925;
int oddrate = 1;
double pass0_frac = 1.0;
int reduce_len = 5;
int reduce_how = 2;
double go_extra = 3.5;
double reduce_margin = 0.06;
int do_reduce = 1;
int pass_threshold = 1;
int strength_how = 1;
int coarse_strength_how = 6;
double reduce_shoulder = -1;
double coarse_all = -1;
int second_count = 2;
int soft_phase_win = 2;
double sub_wrap_thresh = 1.5;
double phase_wrap_thresh = 1.5;
int fancy_subtract = 0;

//
// return a Hamming window of length n.
//
std::vector<double>
hamming(int n)
{
  std::vector<double> h(n);
  for(int k = 0; k < n; k++){
    h[k] = 0.54 - 0.46 * cos(2 * M_PI * k / (n - 1.0));
  }
  return h;
}

//
// blackman window
//
std::vector<double>
blackman(int n)
{
  std::vector<double> h(n);
  for(int k = 0; k < n; k++){
    h[k] = 0.42 - 0.5 * cos(2 * M_PI * k / n) + 0.08*cos(4 * M_PI * k / n);
  }
  return h;
}

//
// symmetric blackman window
//
std::vector<double>
sym_blackman(int n)
{
  std::vector<double> h(n);
  for(int k = 0; k < (n/2)+1; k++){
    h[k] = 0.42 - 0.5 * cos(2 * M_PI * k / n) + 0.08*cos(4 * M_PI * k / n);
  }
  for(int k = n-1; k >= (n/2)+1; --k){
    h[k] = h[(n-1)-k];
  }
  return h;
}

//
// blackman-harris window
//
std::vector<double>
blackmanharris(int n)
{
  double a0 = 0.35875;
  double a1 = 0.48829;
  double a2 = 0.14128;
  double a3 = 0.01168;
  std::vector<double> h(n);
  for(int k = 0; k < n; k++){
    // symmetric
    h[k] =
      a0 
      - a1 * cos(2 * M_PI * k / (n-1))
      + a2 * cos(4 * M_PI * k / (n-1))
      - a3 * cos(6 * M_PI * k / (n-1));
    // periodic
    //h[k] =
    //  a0 
    //  - a1 * cos(2 * M_PI * k / n)
    //  + a2 * cos(4 * M_PI * k / n)
    //  - a3 * cos(6 * M_PI * k / n);
  }
  return h;
}

//
// check the FT8 CRC-14
//
int
check_crc(const int a91[91])
{
  int aa[91];
  int non_zero = 0;
  for(int i = 0; i < 91; i++){
    if(i < 77){
      aa[i] = a91[i];
    } else {
      aa[i] = 0;
    }
    if(aa[i])
      non_zero++;
  }
  int out1[14];

  extern void ft8_crc(int msg1[], int msglen, int out[14]);

  // don't bother with all-zero messages.
  if(non_zero == 0)
    return 0;
  
  // why 82? why not 77?
  ft8_crc(aa, 82, out1);

  for(int i = 0; i < 14; i++){
    if(out1[i] != a91[91-14+i]){
      return 0;
    }
  }
  return 1;
}

//
// manage statistics for soft decoding, to help
// decide how likely each symbol is to be correct,
// to drive LDPC decoding.
//
class Stats {
public:
  std::vector<double> a_;
  bool finalized_;
  double mean_; // cached
  double stddev_; // cached
  
public:
  Stats() : finalized_(false) { }

  void add(double x) {
    a_.push_back(x);
    finalized_ = false;
  }

  void finalize() {
    finalized_ = true;
    
    int n = a_.size();
    mean_ = 0.0;
    for(int i = 0; i < n; i++){
      mean_ += a_[i];
    }
    mean_ /= n;

    double var = 0;
    for(int i = 0; i < n; i++){
      double y = a_[i] - mean_;
      var += y * y;
    }
    var /= n;
    stddev_ = sqrt(var);

    // prepare for binary search to find where values lie
    // in the distribution.
    std::sort(a_.begin(), a_.end());
  }

  double mean() {
    if(!finalized_)
      finalize();
    return mean_;
  }

  double stddev() {
    if(!finalized_)
      finalize();
    return stddev_;
  }

  // fraction of distribution that's less than x.
  // assumes normal distribution.
  double gaussian_problt(double x) {
    double SDs = (x - mean()) / stddev();
    double frac = 0.5 * (1.0 + erf(SDs / sqrt(2.0)));
    return frac;
  }

  // look into the actual distribution.
  double problt(double x, int how) {
    if(!finalized_)
      finalize();

    if(how == 0){
      return gaussian_problt(x);
    }

    // binary search.
    auto it = std::lower_bound(a_.begin(), a_.end(), x);
    int i = it - a_.begin();
    int n = a_.size();

    if(how == 1){
      // index into the distribution.
      // works poorly for values that are off the ends
      // of the distribution, since those are all
      // mapped to 0.0 or 1.0, regardless of magnitude.
      return i / (double) n;
    }

    if(how == 2){
      // use a kind of logistic regression for
      // values near the edges of the distribution.
      if(i < log_tail * n){
        double x0 = a_[(int)(log_tail * n)];
        double y = 1.0 / (1.0 + exp(-log_rate*(x-x0)));
        // y is 0..0.5
        y /= 5;
        return y;
      } else if(i > (1-log_tail) * n){
        double x0 = a_[(int)((1-log_tail) * n)];
        double y = 1.0 / (1.0 + exp(-log_rate*(x-x0)));
        // y is 0.5..1
        // we want (1-log_tail)..1
        y -= 0.5;
        y *= 2;
        y *= log_tail;
        y += (1-log_tail);
        return y;
      } else {
        return i / (double) n;
      }
    }

    if(how == 3){
      // gaussian for values near the edge of the distribution.
      if(i < log_tail * n){
        return gaussian_problt(x);
      } else if(i > (1-log_tail) * n){
        return gaussian_problt(x);
      } else {
        return i / (double) n;
      }
    }

    if(how == 4){
      // gaussian for values outside the distribution.
      if(x < a_[0] || x > a_.back()){
        return gaussian_problt(x);
      } else {
        return i / (double) n;
      }
    }

    assert(0);
  }
};

// a-priori probability of each of the 174 LDPC codeword
// bits being one. measured from reconstructed correct
// codewords, into ft8bits, then python bprob.py.
// from ft8-n4
double apriori174[] = {
  0.47, 0.32, 0.29, 0.37, 0.52, 0.36, 0.40, 0.42, 0.42, 0.53, 0.44,
  0.44, 0.39, 0.46, 0.39, 0.38, 0.42, 0.43, 0.45, 0.51, 0.42, 0.48,
  0.31, 0.45, 0.47, 0.53, 0.59, 0.41, 0.03, 0.50, 0.30, 0.26, 0.40,
  0.65, 0.34, 0.49, 0.46, 0.49, 0.69, 0.40, 0.45, 0.45, 0.60, 0.46,
  0.43, 0.49, 0.56, 0.45, 0.55, 0.51, 0.46, 0.37, 0.55, 0.52, 0.56,
  0.55, 0.50, 0.01, 0.19, 0.70, 0.88, 0.75, 0.75, 0.74, 0.73, 0.18,
  0.71, 0.35, 0.60, 0.58, 0.36, 0.60, 0.38, 0.50, 0.02, 0.01, 0.98,
  0.48, 0.49, 0.54, 0.50, 0.49, 0.53, 0.50, 0.49, 0.49, 0.51, 0.51,
  0.51, 0.47, 0.50, 0.53, 0.51, 0.46, 0.51, 0.51, 0.48, 0.51, 0.52,
  0.50, 0.52, 0.51, 0.50, 0.49, 0.53, 0.52, 0.50, 0.46, 0.47, 0.48,
  0.52, 0.50, 0.49, 0.51, 0.49, 0.49, 0.50, 0.50, 0.50, 0.50, 0.51,
  0.50, 0.49, 0.49, 0.55, 0.49, 0.51, 0.48, 0.55, 0.49, 0.48, 0.50,
  0.51, 0.50, 0.51, 0.50, 0.51, 0.53, 0.49, 0.54, 0.50, 0.48, 0.49,
  0.46, 0.51, 0.51, 0.52, 0.49, 0.51, 0.49, 0.51, 0.50, 0.49, 0.50,
  0.50, 0.47, 0.49, 0.52, 0.49, 0.51, 0.49, 0.48, 0.52, 0.48, 0.49,
  0.47, 0.50, 0.48, 0.50, 0.49, 0.51, 0.51, 0.51, 0.49,
};

class FT8 {
public:
  std::thread *th_;

  double min_hz_;
  double max_hz_;
  std::vector<double> samples_;  // input to each pass
  std::vector<double> nsamples_; // subtract from here

  int start_; // sample number of 0.5 seconds into samples[]
  int rate_;  // samples/second
  double deadline_; // start time + budget
  double final_deadline_; // keep going this long if no decodes
  std::vector<int> hints1_;
  std::vector<int> hints2_;
  int pass_;
  double down_hz_;

  static std::mutex cb_mu_;
  cb_t cb_; // call-back into Python

  std::mutex hack_mu_;
  int hack_size_;
  int hack_off_;
  int hack_len_;
  double hack_0_;
  double hack_1_;
  const double *hack_data_;
  std::vector<std::complex<double>> hack_bins_;

  FT8(const std::vector<double> &samples,
      double min_hz,
      double max_hz,
      int start, int rate,
      int hints1[], int hints2[], double deadline,
      double final_deadline, cb_t cb)
  {
    samples_ = samples;
    min_hz_ = min_hz;
    max_hz_ = max_hz;
    
    start_ = start;
    rate_ = rate;
    deadline_ = deadline;
    final_deadline_ = final_deadline;
    cb_ = cb;
    down_hz_ = 0;

    for(int i = 0; hints1[i]; i++){
      hints1_.push_back(hints1[i]);
    }
    for(int i = 0; hints2[i]; i++){
      hints2_.push_back(hints2[i]);
    }

    hack_size_ = -1;
    hack_data_ = 0;
    hack_off_ = -1;
    hack_len_ = -1;
  }

  ~FT8() {
  }


// strength of costas block of signal with tone 0 at bi0,
// and symbol zero at si0.
double
one_coarse_strength(const ffts_t &bins, int bi0, int si0)
{
  int costas[] = { 3, 1, 4, 0, 6, 5, 2 };

  assert(si0 >= 0 && si0+72+8 <= bins.size());
  assert(bi0 >= 0 && bi0 + 8 <= bins[0].size());

  double sig = 0.0;
  double noise = 0.0;

  if(coarse_all >= 0){
    for(int si = 0; si < 79; si++){
      double mx;
      int mxi = -1;
      double sum = 0;
      for(int i = 0; i < 8; i++){
        double x = std::abs(bins[si0+si][bi0+i]);
        sum += x;
        if(mxi < 0 || x > mx){
          mxi = i;
          mx = x;
        }
      }
      if(si >= 0 && si < 7){
        double x = std::abs(bins[si0+si][bi0+costas[si-0]]);
        sig += x;
        noise += sum - x;
      } else if(si >= 36 && si < 36+7){
        double x = std::abs(bins[si0+si][bi0+costas[si-36]]);
        sig += x;
        noise += sum - x;
      } else if(si >= 72 && si < 72+7){
        double x = std::abs(bins[si0+si][bi0+costas[si-72]]);
        sig += x;
        noise += sum - x;
      } else {
        sig += coarse_all * mx;
        noise += coarse_all * (sum - mx);
      }
    }
  } else {
    // just costas symbols
    for(int si = 0; si < 7; si++){
      for(int bi = 0; bi < 8; bi++){
        double x = 0;
        x += std::abs(bins[si0+si][bi0+bi]);
        x += std::abs(bins[si0+36+si][bi0+bi]);
        x += std::abs(bins[si0+72+si][bi0+bi]);
        if(bi == costas[si]){
          sig += x;
        } else {
          noise += x;
        }
      }
    }
  }

  if(coarse_strength_how == 0){
     return sig - noise;
  } else if(coarse_strength_how == 1){
    return sig - noise/7;
  } else if(coarse_strength_how == 2){
    return sig / (noise/7);
  } else if(coarse_strength_how == 3){
    return sig / (sig + (noise/7));
  } else if(coarse_strength_how == 4){
    return sig;
  } else if(coarse_strength_how == 5){
    return sig / (sig + noise);
  } else if(coarse_strength_how == 6){
    return sig / noise;
  } else {
    assert(0);
  }
}

// return symbol length in samples at the given rate.
// insist on integer symbol lengths so that we can
// use whole FFT bins.
int
blocksize(int rate)
{
  // FT8 symbol length is 1920 at 12000 samples/second.
  int xblock = 1920 / (12000.0 / rate);
  assert(xblock == (int) xblock);
  int block = xblock;
  return block;
}

class Strength {
public:
  double hz_;
  int off_;
  double strength_; // higher is better
};

//
// look for potential signals by searching FFT bins for Costas symbol
// blocks. returns a vector of candidate positions.
//
std::vector<Strength>
coarse(const ffts_t &bins, int si0, int si1)
{
  int block = blocksize(rate_);
  int nbins = bins[0].size();
  double bin_hz = rate_ / (double) block;
  int min_bin = min_hz_ / bin_hz;
  int max_bin = max_hz_ / bin_hz;

  std::vector<Strength> strengths;
  
  for(int bi = min_bin; bi < max_bin && bi+8 <= nbins; bi++){
    std::vector<Strength> sv;
    for(int si = si0; si < si1 && si + 79 < bins.size(); si++){
      double s = one_coarse_strength(bins, bi, si);
      Strength st;
      st.strength_ = s;
      st.hz_ = bi * 6.25;
      st.off_ = si * block;
      sv.push_back(st);
    }
    if(sv.size() < 1)
      break;

    // save best ncoarse offsets, but require that they be separated
    // by at least one symbol time.

    std::sort(sv.begin(), sv.end(),
              [](const Strength &a, const Strength &b) -> bool
              { return a.strength_ > b.strength_; } );

    strengths.push_back(sv[0]);

    int nn = 1;
    for(int i = 1; nn < ncoarse && i < sv.size(); i++){
      if(std::abs(sv[i].off_ - sv[0].off_) > ncoarse_blocks*block){
        strengths.push_back(sv[i]);
        nn++;
      }
    }
  }

  return strengths;
}

//
// reduce the sample rate from arate to brate.
// center hz0..hz1 in the new nyquist range.
// but first filter to that range.
// sets delta_hz to hz moved down.
//
std::vector<double>
reduce_rate(std::vector<double> a, double hz0, double hz1,
            int arate, int brate,
            double &delta_hz)
{
  assert(brate < arate);
  assert(hz1 - hz0 <= brate / 2);

  // the pass band is hz0..hz1
  // stop bands are 0..hz00 and hz11..nyquist.
  double hz00, hz11;

  if(reduce_shoulder > 0){
    hz00 = hz0 - reduce_shoulder;
    hz11 = hz1 + reduce_shoulder;
  } else {
    double mid = (hz0 + hz1) / 2;
    hz00 = mid - (brate / 4.0);
    hz11 = mid + (brate / 4.0);
  }

  int len = a.size();
  std::vector<std::complex<double>> bins = one_fft(a, 0, len, "reduce_rate");
  int nbins = bins.size();
  double bin_hz = arate / (double) len;

  if(reduce_how == 2){
    // band-pass filter the FFT output.
    bins = fbandpass(bins, bin_hz,
                     hz00,
                     hz0,
                     hz1,
                     hz11);
  }
  
  if(reduce_how == 3){
    for(int i = 0; i < nbins; i++){
      if(i < (hz0 / bin_hz)){
        bins[i] = 0;
      } else if(i > (hz1 / bin_hz)){
        bins[i] = 0;
      }
    }
  }

  // shift down.
  int omid = ((hz0 + hz1) / 2) / bin_hz;
  int nmid = (brate / 4) / bin_hz;

  int delta = omid - nmid; // amount to move down
  assert(delta < nbins);
  if(delta > 0){
    for(int i = 0; i < nbins-delta; i++){
      bins[i] = bins[i + delta];
    }
    for(int i = nbins-delta; i < nbins; i++){
      bins[i] = 0;
    }
  }

  // use ifft to reduce the rate.
  int blen = round(a.size() * (brate / (double) arate));
  std::vector<std::complex<double>> bbins(blen / 2 + 1);
  for(int i = 0; i < bbins.size(); i++)
    bbins[i] = bins[i];
  std::vector<double> vvv = one_ifft(bbins, "reduce_rate");

  delta_hz = delta * bin_hz;

  return vvv;
}

void
go()
{
#if 0
  fprintf(stderr, "go: %.0f .. %.0f, %.0f, rate=%d\n",
          min_hz_, max_hz_, max_hz_ - min_hz_, rate_);
#endif

  // can we reduce the sample rate?
  int rates[] = { 500, 1000, 1500, 2000, 3000, 4000, 6000, -1 };
  int nrate = -1;
#if 0
  for(int ratei = 0; rates[ratei] > 0; ratei++){
    int xrate = rates[ratei];
#else
  //for(int xrate = 25; xrate < rate_; xrate += 25){
  for(int xrate = 100; xrate < rate_; xrate += 100){
#endif
    if(xrate < rate_ && (oddrate || (rate_ % xrate) == 0)){
      if(((max_hz_ - min_hz_) + 50 + 2*go_extra) < nyquist * (xrate / 2)){
        nrate = xrate;
        break;
      }
    }
  }

  if(do_reduce && nrate > 0 && nrate < rate_){
    // filter and reduce the sample rate from rate_ to nrate.

    double delta_hz; // how much it moved down
    samples_ = reduce_rate(samples_,
                           min_hz_-3.1-go_extra,
                           max_hz_+50-3.1+go_extra,
                           rate_, nrate, delta_hz);

    if(0){
      fprintf(stderr, "%.0f..%.0f, range %.0f, rate %d -> %d, delta hz %.0f\n",
              min_hz_, max_hz_,
              max_hz_ - min_hz_,
              rate_, nrate, delta_hz);
    }

    if(delta_hz > 0){
      down_hz_ = delta_hz; // to adjust hz for Python.
      min_hz_ -= down_hz_;
      max_hz_ -= down_hz_;
    }
    assert(max_hz_ + 50 < nrate / 2);
    assert(min_hz_ >= 0);

    double ratio = nrate / (double) rate_;
    rate_ = nrate;
    start_ = round(start_ * ratio);
  }
  
  int block = blocksize(rate_);

  double bin_hz = rate_ / (double) block;

  // start_ is the sample number of 0.5 seconds, the nominal start time.

#if 0
  // make sure there's at least tminus*rate_ samples before start_.
  // XXX this breaks reported offsets.
  if(start_ < tminus*rate_){
    int need = tminus*rate_ - start_;
    fprintf(stderr, "*** prepending %d\n", need);
    std::vector<double> v(need);
    for(int i = 0; i < need; i++)
      v[i] = 0;
    samples_.insert(samples_.begin(), v.begin(), v.end());
    start_ += need;
  }
#endif

  // make sure there's at least tplus*rate_ samples after the end.
  if(start_ + tplus*rate_ + 79 * block > samples_.size()){
    int need = start_ + tplus*rate_ + 79 * block - samples_.size();

    // round up to a whole second, to ease fft plan caching.
    if((need % rate_) != 0)
      need += rate_ - (need % rate_);

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, samples_.size()-1);
    auto rnd = std::bind(distribution, generator);

    std::vector<double> v(need);
    for(int i = 0; i < need; i++){
      //v[i] = 0;
      v[i] = samples_[rnd()];
    }
    samples_.insert(samples_.end(), v.begin(), v.end());
  }
  
  int si0 = (start_ - tminus*rate_) / block;
  if(si0 < 0)
    si0 = 0;
  int si1 = (start_ + tplus*rate_) / block;
  
  // a copy from which to subtract.
  nsamples_ = samples_;

  for(pass_ = 0; pass_ < npasses; pass_++){
    double total_remaining = deadline_ - now();
    double remaining = total_remaining / (npasses - pass_);
    if(pass_ == 0){
      remaining *= pass0_frac;
    }
    double deadline = now() + remaining;
    
    int new_decodes = 0;
    samples_ = nsamples_;

    std::vector<Strength> order;

    //
    // search coarsely for Costas blocks.
    // in fractions of bins in off and hz.
    //

    for(int hz_frac_i = 0; hz_frac_i < coarse_hz_fracs; hz_frac_i++){
      // shift down by hz_frac
      double hz_frac = hz_frac_i * (6.25 / coarse_hz_fracs);
      std::vector<double> samples1;
      if(hz_frac_i == 0){
        samples1 = samples_;
      } else {
        samples1 = fft_shift(samples_, 0, samples_.size(),
                             rate_, hz_frac);
      }
      
      for(int off_frac_i = 0; off_frac_i < coarse_off_fracs; off_frac_i++){
        int off_frac = off_frac_i * (block / coarse_off_fracs);
        ffts_t bins = ffts(samples1, off_frac, block, "go");
        std::vector<Strength> oo = coarse(bins, si0, si1);
        for(int i = 0; i < oo.size(); i++){
          oo[i].hz_ += hz_frac;
          oo[i].off_ += off_frac;
        }
        order.insert(order.end(), oo.begin(), oo.end());
      }
    }

    //
    // sort strongest-first.
    //
    std::sort(order.begin(), order.end(),
              [](const Strength &a, const Strength &b) -> bool
              { return a.strength_ > b.strength_; } );
    
    char already[2000]; // XXX
    for(int i = 0; i < sizeof(already)/sizeof(already[0]); i++)
      already[i] = 0;
    
    for(int ii = 0; ii < order.size(); ii++){
      double tt = now();
      if(ii > 0 &&
         tt > deadline &&
         (tt > deadline_ || new_decodes >= pass_threshold) &&
         (pass_ < npasses-1 || tt > final_deadline_)){
        break;
      }

      double hz = order[ii].hz_;
      if(already[(int)round(hz / already_hz)])
        continue;
      int off = order[ii].off_;
      int ret = one(samples_, hz, off);
      if(ret){
        if(ret == 2){
          new_decodes++;
        }
        already[(int)round(hz / already_hz)] = 1;
      }
    }
  }
}

//
// what's the strength of the Costas sync blocks of
// the signal starting at hz and off?
//
double
one_strength(const std::vector<double> &samples200, double hz, int off)
{
  int bin0 = round(hz / 6.25);

  int costas[] = { 3, 1, 4, 0, 6, 5, 2 };

  double sig = 0;
  double noise = 0;

  for(int si = 0; si < 7; si++){
    auto fft1 = one_fft(samples200, off+(si+0)*32, 32, "one_strength");
    auto fft2 = one_fft(samples200, off+(si+36)*32, 32, "one_strength");
    auto fft3 = one_fft(samples200, off+(si+72)*32, 32, "one_strength");
    for(int bi = 0; bi < 8; bi++){
      double x = 0;
      x += std::abs(fft1[bin0+bi]);
      x += std::abs(fft2[bin0+bi]);
      x += std::abs(fft3[bin0+bi]);
      if(bi == costas[si]){
        sig += x;
      } else {
        noise += x;;
      }
    }
  }
  
  if(strength_how == 0){
     return sig - noise;
  } else if(strength_how == 1){
    return sig - noise/7;
  } else if(strength_how == 2){
    return sig / (noise/7);
  } else if(strength_how == 3){
    return sig / (sig + (noise/7));
  } else if(strength_how == 4){
    return sig;
  } else if(strength_how == 5){
    return sig / (sig + noise);
  } else if(strength_how == 6){
    return sig / noise;
  } else {
    assert(0);
  }
}

//
// given a complete known signal's symbols in syms,
// how strong is it? used to look for the best
// offset and frequency at which to subtract a
// decoded signal.
//
double
one_strength_known(const std::vector<double> &samples200,
                   const std::vector<int> syms,
                   double hz, int off)
{
  assert(syms.size() == 79);
  
  int bin0 = round(hz / 6.25);

  double sig = 0;
  double noise = 0;
  for(int si = 0; si < 79; si++){
    auto fft1 = one_fft(samples200, off+si*32, 32, "one_strength_known");
    for(int bi = 0; bi < 8; bi++){
      double x = std::abs(fft1[bin0+bi]);
      if(bi == syms[si]){
        sig += x;
      } else {
        noise += x;
      }
    }
  }
  if(strength_how == 0){
     return sig - noise;
  } else if(strength_how == 1){
    return sig - noise/7;
  } else if(strength_how == 2){
    return sig / (noise/7);
  } else if(strength_how == 3){
    return sig / (sig + (noise/7));
  } else if(strength_how == 4){
    return sig;
  } else if(strength_how == 5){
    return sig / (sig + noise);
  } else if(strength_how == 6){
    return sig / noise;
  } else {
    assert(0);
  }
}

int
search_time_fine(const std::vector<double> &samples200,
                 int off0, int offN,
                 double hz,
                 int gran,
                 double &str)
{
  if(off0 < 0)
    off0 = 0;

  //
  // shift in frequency to put hz at 25.
  // only shift the samples we need, both for speed,
  // and try to always shift down the same number of samples
  // to make it easier to cache fftw plans.
  //
  int len = (offN - off0) + 79*32 + 32;
  if(off0 + len > samples200.size()){
    // len = samples200.size() - off0;
    // don't provoke random-length FFTs.
    return -1;
  }
  std::vector<double> downsamples200 = shift200(samples200, off0, len, hz);

  int best_off = -1;
  double best_sum = 0.0;

  for(int g = 0; g <= (offN-off0) && g + 79*32 <= len; g += gran){
    double sum = one_strength(downsamples200, 25, g);
    if(sum > best_sum || best_off == -1){
      best_off = g;
      best_sum = sum;
    }
  }

  str = best_sum;
  assert(best_off >= 0);
  return off0 + best_off;
}

int
search_time_fine_known(const std::vector<double> &samples200,
                       const std::vector<int> &syms,
                       int off0, int offN,
                       double hz,
                       int freq_factor, int time_factor,
                       double &str)
{
  if(off0 < 0)
    off0 = 0;
  
  // put hz at 25.
  std::vector<double> downsamples200 = shift200(samples200, 0, samples200.size(), hz);

  assert(time_factor > 0);
  int gran = 32 / time_factor;

  int best_off = -1;
  double best_sum = 0.0;
  int g = off0;
  if(g < 0)
    g = 0;

  for( ; g <= offN && g + 79*32 <= downsamples200.size(); g += gran){
    double sum = one_strength_known(downsamples200, syms, 25, g);
    if(sum > best_sum || best_off == -1){
      best_off = g;
      best_sum = sum;
    }
  }

  if(best_off < 0)
    return -1;

  str = best_sum;
  return best_off;
}

//
// search for costas blocks in an MxN time/frequency grid.
// hz0 +/- hz_win in hz_inc increments. hz0 should be near 25.
// off0 +/- off_win in off_inc incremenents.
//
std::vector<Strength>
search_both(const std::vector<double> &samples200,
            double hz0, int hz_fracs, double hz_win,
            int off0, int off_fracs, int off_win)
{
  assert(hz0 >= 25 - 6.25/2 && hz0 <= 25 + 6.25/2);

  std::vector<Strength> strengths;
  
  double hz_inc = 2 * hz_win / hz_fracs;
  int off_inc = round(2 * off_win / (double) off_fracs);
  

  for(double hz = hz0 - hz_win; hz <= hz0 + hz_win + 0.01; hz += hz_inc){
    double str = 0;
    int off = search_time_fine(samples200, off0 - off_win, off0 + off_win, hz,
                               off_inc, str);
    if(off >= 0){
      Strength st;
      st.hz_ = hz;
      st.off_ = off;
      st.strength_ = str;
      strengths.push_back(st);
    }
  }

  return strengths;
}

void
search_both_known(const std::vector<double> &samples200,
                  const std::vector<int> &syms,
                  double hz0, double hz_inc, double hz_win,
                  int off0, int off_inc, int off_win,
                  double &hz_out, int &off_out, double &strength_out)
{
  assert(hz0 >= 25 - 6.25 && hz0 <= 25 + 6.25);
  
  int got_best = 0;
  double best_hz = 0;
  int best_off = 0;
  double best_strength = 0;

  for(double hz = hz0 - hz_win; hz <= hz0 + hz_win; hz += hz_inc){
    double strength = 0;
    int time_factor = 32 / off_inc;
    int off = search_time_fine_known(samples200, syms,
                                     off0 - off_win, off0 + off_win, hz,
                                     0, time_factor, strength);
    if(off >= 0 && (got_best == 0 || strength > best_strength)){
      got_best = 1;
      best_hz = hz;
      best_off = off;
      best_strength = strength;
    }
  }

  if(got_best){
    hz_out = best_hz;
    off_out = best_off;
    strength_out = best_strength;
  }
}

//
// shift frequency by shifting the bins of one giant FFT.
// so no problem with phase mismatch &c at block boundaries.
// surprisingly fast at 200 samples/second.
// shifts *down* by hz.
//
std::vector<double>
fft_shift(const std::vector<double> &samples, int off, int len,
          int rate, double hz)
{
  // horrible hack to avoid repeated FFTs on the same input.
  hack_mu_.lock();
  std::vector<std::complex<double>> bins;
  if(samples.size() == hack_size_ && samples.data() == hack_data_ &&
     off == hack_off_ && len == hack_len_ &&
     samples[0] == hack_0_ && samples[1] == hack_1_){
    bins = hack_bins_;
  } else {
    bins = one_fft(samples, off, len, "fft_shift");
    hack_bins_ = bins;
    hack_size_ = samples.size();
    hack_off_ = off;
    hack_len_ = len;
    hack_0_ = samples[0];
    hack_1_ = samples[1];
    hack_data_ = samples.data();
  }
  hack_mu_.unlock();
  
  int nbins = bins.size();

  double bin_hz = rate / (double) len;
  int down = round(hz / bin_hz);
  std::vector<std::complex<double>> bins1(nbins);
  for(int i = 0; i < nbins; i++){
    int j = i + down;
    if(j >= 0 && j < nbins){
      bins1[i] = bins[j];
    } else {
      bins1[i] = 0;
    }
  }
  std::vector<double> out = one_ifft(bins1, "fft_shift");
  return out;
}

// shift the frequency by a fraction of 6.25,
// to center hz on bin 4 (25 hz).
std::vector<double>
shift200(const std::vector<double> &samples200, int off, int len, double hz)
{
  if(std::abs(hz - 25) < 0.001 && off == 0 && len == samples200.size()){
    return samples200;
  } else {
    return fft_shift(samples200, off, len, 200, hz - 25.0);
  }
  // return hilbert_shift(samples200, hz - 25.0, hz - 25.0, 200);
}

// returns a mini-FFT of 79 8-tone symbols.
ffts_t
extract(const std::vector<double> &samples200, double hz, int off)
{
  // put hz in the middle of bin 4, at 25.0.
  std::vector<double> downsamples200 = shift200(samples200, 0, samples200.size(), hz);

  ffts_t bins3 = ffts(downsamples200, off, 32, "extract");

  ffts_t m79(79);
  for(int si = 0; si < 79; si++){
    m79[si].resize(8);
    if(si < (int)bins3.size()){
      for(int bi = 0; bi < 8; bi++){
        auto x = bins3[si][4+bi];
        m79[si][bi] = x;
      }
    } else {
      for(int bi = 0; bi < 8; bi++){
        m79[si][bi] = 0;
      }
    }
  }

  return m79;
}

//
// m79 is a 79x8 array of complex.
//
ffts_t
un_gray_code_c(const ffts_t &m79)
{
  ffts_t m79a(79);

  int map[] = { 0, 1, 3, 2, 6, 4, 5, 7 };
  for(int si = 0; si < 79; si++){
    m79a[si].resize(8);
    for(int bi = 0; bi < 8; bi++){
      m79a[si][map[bi]] = m79[si][bi];
    }
  }

  return m79a;
}

//
// m79 is a 79x8 array of double.
//
std::vector< std::vector<double> > 
un_gray_code_r(const std::vector<std::vector<double>> &m79)
{
  std::vector< std::vector<double> > m79a(79);

  int map[] = { 0, 1, 3, 2, 6, 4, 5, 7 };
  for(int si = 0; si < 79; si++){
    m79a[si].resize(8);
    for(int bi = 0; bi < 8; bi++){
      m79a[si][map[bi]] = m79[si][bi];
    }
  }

  return m79a;
}

//
// turn 79 symbol numbers into 174 bits.
// strip out the three Costas sync blocks,
// leaving 58 symbol numbers.
// each represents three bits.
// (all post-un-gray-code).
//
std::vector<int>
extract_bits(const std::vector<int> &syms)
{
  assert(syms.size() == 79);

  std::vector<int> bits;
  for(int si = 0; si < 79; si++){
    if(si < 7 || (si >= 36 && si < 36+7) || si >= 72){
      // costas -- skip
    } else {
      bits.push_back((syms[si] & 4) != 0);
      bits.push_back((syms[si] & 2) != 0);
      bits.push_back((syms[si] & 1) != 0);
    }
  }

  return bits;
}

//
// normalize levels by windowed median.
// this helps, but why?
//
std::vector< std::vector<double> >
convert_to_snr(const std::vector< std::vector<double> > &m79, int how, int win)
{
  if(how < 0 || win < 0)
    return m79;

  // for each symbol time, what's its "noise" level?
  //
  std::vector<double> mm(79);
  for(int si = 0; si < 79; si++){
    std::vector<double> v(8);
    double sum = 0.0;
    for(int bi = 0; bi < 8; bi++){
      double x = m79[si][bi];
      v[bi] = x;
      sum += x;
    }
    if(how != 1)
      std::sort(v.begin(), v.end());
    if(how == 0){
      // median
      mm[si] = (v[3] + v[4]) / 2;
    } else if(how == 1){
      mm[si] = sum / 8;
    } else if(how == 2){
      // all but strongest tone.
      mm[si] = (v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6]) / 7;
    } else if(how == 3){
      mm[si] = v[0]; // weakest tone
    } else if(how == 4){
      mm[si] = v[7]; // strongest tone
    } else if(how == 5){
      mm[si] = v[6]; // second-strongest tone
    } else {
      mm[si] = 1.0;
    }
  }

  // we're going to take a windowed average.
  std::vector<double> winwin;
  if(win > 0){
    winwin = blackman(2*win+1);
  } else {
    winwin.push_back(1.0);
  }

  std::vector<std::vector<double>> n79(79);
    
  for(int si = 0; si < 79; si++){
    double sum = 0;
    for(int dd = si - win; dd <= si + win; dd++){
      int wi = dd - (si - win);
      if(dd >= 0 && dd < 79){
        sum += mm[dd] * winwin[wi];
      } else if(dd < 0){
        sum += mm[0] * winwin[wi];
      } else {
        sum += mm[78] * winwin[wi];
      }
    }
    n79[si].resize(8);
    for(int bi = 0; bi < 8; bi++){
      n79[si][bi] = m79[si][bi] / sum;
    }
  }

  return n79;
}

//
// statistics to decide soft probabilities,
// to drive LDPC decoder.
// distribution of strongest tones, and
// distribution of noise.
// multiple ranges in case things change over time.
// nranges is soft_ranges.
//
void
make_stats(const std::vector<std::vector<double>> &m79,
           std::vector<Stats> &noises,
           std::vector<Stats> &bests,
           int nranges,
           int x_best_in_noise)
{
  noises.resize(nranges);
  bests.resize(nranges);

  int costas[] = { 3, 1, 4, 0, 6, 5, 2 };

  for(int range = 0; range < nranges; range++){
    int si0 = range * (79 / nranges);
    int si1;
    if(range == nranges - 1)
      si1 = 79;
    else
      si1 = (range + 1) * (79 / nranges);
    for(int si = si0; si < si1 && si < 79; si++){
      if(si < 7 || (si >= 36 && si < 36 + 7) || si >= 72){
        // Costas.
        int ci;
        if(si >= 72) ci = si - 72;
        else if(si >= 36) ci = si - 36;
        else ci = si;
        for(int bi = 0; bi < 8; bi++){
          double x = m79[si][bi];
          if(bi == costas[ci]){
            bests[range].add(x);
            if(x_best_in_noise)
              noises[range].add(x); // include best in noises too
          } else {
            noises[range].add(x);
          }
        }
      } else {
        std::vector<double> v(8);
        for(int bi = 0; bi < 8; bi++){
          v[bi] = m79[si][bi];
        }
        std::sort(v.begin(), v.end());
        for(int i = 0; i < 7; i++){
          noises[range].add(v[i]);
        }
        bests[range].add(v[7]);
        if(x_best_in_noise)
          noises[range].add(v[7]); // include best in noises too
      }
    }
  }
}

//
// convert 79x8 complex FFT bins to magnitudes.
//
// exploits local phase coherence by decreasing magnitudes of bins
// whose phase is far from the phases of nearby strongest tones.
//
// relies on each tone being reasonably well centered in its FFT bin
// (in time and frequency) so that each tone completes an integer
// number of cycles and thus preserves phase from one symbol to the
// next.
//
std::vector<std::vector<double>>
soft_c2m(const ffts_t &c79)
{
  std::vector< std::vector<double> > m79(79);
  std::vector<double> raw_phases(79); // of strongest tone in each symbol time
  for(int si = 0; si < 79; si++){
    m79[si].resize(8);
    int mxi = -1;
    double mx;
    double mx_phase;
    for(int bi = 0; bi < 8; bi++){
      double x = std::abs(c79[si][bi]);
      m79[si][bi] = x;
      if(mxi < 0 || x > mx){
        mxi = bi;
        mx = x;
        mx_phase = std::arg(c79[si][bi]); // -pi .. pi
      }
    }
    raw_phases[si] = mx_phase;
  }

  if(soft_phase_win <= 0)
    return m79;

  // phase around each symbol.
  std::vector<double> phases(79);

  for(int si = 0; si < 79; si++){
    // median of nearby phases
    std::vector<double> v;
    int nhi = 0;
    int nlo = 0;
    for(int si1 = si - soft_phase_win; si1 <= si + soft_phase_win; si1++){
      if(si1 >= 0 && si1 < 79){
        double x = raw_phases[si1];
        v.push_back(x);
#if 0
        if(x < -2)
          nlo++;
        if(x > 2)
          nhi++;
#else
        if(x < 0)
          nlo++;
        else
          nhi++;
#endif
      }
    }

    // try to fix -pi / pi wrap-around.
    int n = v.size();
    if(nhi + nlo >= n-1){
      if(nhi > nlo){
        for(int i = 0; i < n; i++){
          if(v[i] < -phase_wrap_thresh){
            v[i] += 2*M_PI;
          }
        }
      } else if(nlo > nhi){
        for(int i = 0; i < n; i++){
          if(v[i] > phase_wrap_thresh){
            v[i] -= 2*M_PI;
          }
        }
      }
    }
       
    std::sort(v.begin(), v.end());
    phases[si] = v[n/2];
  }

  // project each tone against the median phase around that symbol time.
  for(int si = 0; si < 79; si++){
    for(int bi = 0; bi < 8; bi++){
      double mag = std::abs(c79[si][bi]);
      double angle = std::arg(c79[si][bi]);
      double d = angle - phases[si];
      double factor = 0.1;
      if(d < M_PI/2 && d > -M_PI/2){
        factor = cos(d);
      }
      m79[si][bi] = factor * mag;
    }
  }
  
  return m79;
}

//
// c79 is 79x8 complex tones, before un-gray-coding.
//
void
prepare_soft(const ffts_t &c79, double ll174[], int best_off)
{
  // m79 = absolute values of c79.
  // still pre-un-gray-coding so we know which
  // are the correct Costas tones.
  std::vector< std::vector<double> > m79(79);
  m79 = soft_c2m(c79);

  m79 = convert_to_snr(m79, snr_how, snr_win);
    
  // statistics to decide soft probabilities.
  // distribution of strongest tones, and
  // distribution of noise.
  // multiple ranges in case things change over time.
  std::vector<Stats> noises;
  std::vector<Stats> bests;
  make_stats(m79, noises, bests, soft_ranges, best_in_noise);

  m79 = un_gray_code_r(m79);

  int lli = 0;
  for(int i79 = 0; i79 < 79; i79++){
    if(i79 < 7 || (i79 >= 36 && i79 < 36+7) || i79 >= 72){
      // Costas, skip
      continue;
    }

    int range = i79 / (79 / soft_ranges);

    // for each of the three bits, look at the strongest tone
    // that would make it a zero, and the strongest tone that
    // would make it a one. use Bayes to decide which is more
    // likely, comparing each against the distribution of noise
    // and the distribution of strongest tones.
    // most-significant-bit first.

    for(int biti = 0; biti < 3; biti++){
      // tone numbers that make this bit zero or one.
      int zeroi[4];
      int onei[4];
      if(biti == 0){
        // high bit
        zeroi[0] = 0; zeroi[1] = 1; zeroi[2] = 2; zeroi[3] = 3;
        onei[0] = 4; onei[1] = 5; onei[2] = 6; onei[3] = 7;
      }
      if(biti == 1){
        // middle bit
        zeroi[0] = 0; zeroi[1] = 1; zeroi[2] = 4; zeroi[3] = 5;
        onei[0] = 2; onei[1] = 3; onei[2] = 6; onei[3] = 7;
      }
      if(biti == 2){
        // low bit
        zeroi[0] = 0; zeroi[1] = 2; zeroi[2] = 4; zeroi[3] = 6;
        onei[0] = 1; onei[1] = 3; onei[2] = 5; onei[3] = 7;
      }

      // strongest tone that would make this bit be zero.
      int got_best_zero = 0;
      double best_zero = 0;
      for(int i = 0; i < 4; i++){
        double x = m79[i79][zeroi[i]];
        if(got_best_zero == 0 || x > best_zero){
          got_best_zero = 1;
          best_zero = x;
        }
      }

      // strongest tone that would make this bit be one.
      int got_best_one = 0;
      double best_one = 0;
      for(int i = 0; i < 4; i++){
        double x = m79[i79][onei[i]];
        if(got_best_one == 0 || x > best_one){
          got_best_one = 1;
          best_one = x;
        }
      }

      double maxlog = 4.97;
      double ll = 0;

      //
      // Bayes combining rule normalization from:
      // http://cs.wellesley.edu/~anderson/writing/naive-bayes.pdf
      //
      // a = P(zero)P(e0|zero)P(e1|zero)
      // b = P(one)P(e0|one)P(e1|one)
      // p = a / (a + b)
      //
      
      double pzero = 0.5;
      double pone = 0.5;
      if(use_apriori){
        pzero = 1.0 - apriori174[lli];
        pone = apriori174[lli];
      }
      
      // zero
      double a = pzero *
        bests[range].problt(best_zero, problt_how) *
        (1.0 - noises[range].problt(best_one, problt_how));
      
      // one
      double b = pone *
        bests[range].problt(best_one, problt_how) *
        (1.0 - noises[range].problt(best_zero, problt_how));
      
      double p;
      if(a + b == 0){
        p = 0.5;
      } else {
        p = a / (a + b);
      }
      
      if(1 - p == 0.0){
        ll = maxlog;
      } else {
        ll = log(p / (1 - p));
      }
        
      if(ll > maxlog)
        ll = maxlog;
      if(ll < -maxlog)
        ll = -maxlog;
      
      ll174[lli++] = ll;
    }
  }
  assert(lli == 174);
}

//
// given log likelyhood for each bit, try LDPC and OSD decoders.
//
int
decode(const double ll174[], int a174[], int use_osd, std::string &comment)
{
  void ldpc_decode(double llcodeword[], int iters, int plain[], int *ok);
  void ldpc_decode_log(double codeword[], int iters, int plain[], int *ok);

  int plain[174]; // will be 0/1 bits.
  int ldpc_ok = 0;     // 83 will mean success.

  ldpc_decode((double*)ll174, ldpc_iters, plain, &ldpc_ok);

  int ok_thresh = 83; // 83 is perfect
  if(ldpc_ok >= ok_thresh){
    // plain[] is 91 systematic data bits, 83 parity bits.
    for(int i = 0; i < 174; i++){
      a174[i] = plain[i];
    }
    if(check_crc(a174)){
      // success!
      return 1;
    }
  }

  if(use_osd && osd_depth >= 0 && ldpc_ok >= osd_ldpc_thresh){
    extern int osd_decode(double codeword[174], int depth, int out[91], int*);
    extern void ldpc_encode(int plain[91], int codeword[174]);

    int oplain[91];
    int got_depth = -1;
    int osd_ok = osd_decode((double*)ll174, osd_depth, oplain, &got_depth);
    if(osd_ok){
      // reconstruct all 174.
      comment += "OSD-" + std::to_string(got_depth) + "-" + std::to_string(ldpc_ok);
      ldpc_encode(oplain, a174);
      return 1;
    }
  }
  
  return 0;
}

//
// bandpass filter some FFT bins.
// smooth transition from stop-band to pass-band,
// so that it's not a brick-wall filter, so that it
// doesn't ring.
//
std::vector<std::complex<double>>
fbandpass(const std::vector<std::complex<double>> &bins0,
          double bin_hz,
          double low_outer,  // start of transition
          double low_inner,  // start of flat area
          double high_inner, // end of flat area
          double high_outer) // end of transition
{
  // assert(low_outer >= 0);
  assert(low_outer <= low_inner);
  assert(low_inner <= high_inner);
  assert(high_inner < high_outer);
  // assert(high_outer <= bin_hz * bins0.size());
  
  int nbins = bins0.size();
  std::vector<std::complex<double>> bins1(nbins);

  for(int i = 0; i < nbins; i++){
    double ihz = i * bin_hz;
    // cos(x)+flat+cos(x) taper
    double factor;
    if(ihz <= low_outer || ihz >= high_outer){
      factor = 0;
    } else if(ihz >= low_outer && ihz < low_inner){
      // rising shoulder
      double theta = (ihz - low_outer) / (low_inner-low_outer); // 0 .. 1
      theta -= 1; // -1 .. 0
      theta *= 3.14159; // -pi .. 0
      factor = cos(theta); // -1 .. 1
      factor = (factor + 1) / 2; // 0 .. 1
    } else if(ihz > high_inner && ihz <= high_outer){
      // falling shoulder
      double theta =  (high_outer - ihz) / (high_outer-high_inner); // 1 .. 0
      theta = 1.0 - theta; // 0 .. 1
      theta *= 3.14159; // 0 .. pi
      factor = cos(theta); // 1 .. -1
      factor = (factor + 1) / 2; // 1 .. 0
    } else {
      factor = 1.0;
    }
    bins1[i] = bins0[i] * factor;
  }

  return bins1;
}

//
// move hz down to 25, filter+convert to 200 samples/second.
//
// like fft_shift(). one big FFT, move bins down and
// zero out those outside the band, then IFFT,
// then re-sample.
//
// XXX maybe merge w/ fft_shift() / shift200().
//
std::vector<double>
down_v7(const std::vector<double> &samples, double hz)
{
  int len = samples.size();
  std::vector<std::complex<double>> bins = one_fft(samples, 0, len, "down_v7");
  int nbins = bins.size();

  double bin_hz = rate_ / (double) len;
  int down = round((hz - 25) / bin_hz);
  std::vector<std::complex<double>> bins1(nbins);
  for(int i = 0; i < nbins; i++){
    int j = i + down;
    if(j >= 0 && j < nbins){
      bins1[i] = bins[j];
    } else {
      bins1[i] = 0;
    }
  }

  // now filter to fit in 200 samples/second.

  double low_inner = 25.0 - shoulder200_extra;
  double low_outer = low_inner - shoulder200;
  if(low_outer < 0)
    low_outer = 0;
  double high_inner = 75 - 6.25 + shoulder200_extra;
  double high_outer = high_inner + shoulder200;
  if(high_outer > 100)
    high_outer = 100;

  bins1 = fbandpass(bins1, bin_hz,
                    low_outer, low_inner,
                    high_inner, high_outer);

  // convert back to time domain and down-sample to 200 samples/second.
  int blen = round(samples.size() * (200.0 / rate_));
  std::vector<std::complex<double>> bbins(blen / 2 + 1);
  for(int i = 0; i < bbins.size(); i++)
    bbins[i] = bins1[i];
  std::vector<double> out = one_ifft(bbins, "down_v7");
  
  return out;
}

//
// putative start of signal is at hz and symbol si0.
// 
// return 2 if it decodes to a brand-new message.
// return 1 if it decodes but we've already seen it,
//   perhaps in a different pass.
// return 0 if we could not decode.
//
// XXX merge with one_iter().
//
int
one(const std::vector<double> &samples, double hz, int off)
{
  //
  // set up to search for best frequency and time offset.
  //

  //
  // move down to 25 hz and re-sample to 200 samples/second,
  // i.e. 32 samples/symbol.
  //
  std::vector<double> samples200 = down_v7(samples, hz);

  int off200 = (off / (double) rate_) * 200;

  int ret = one_iter(samples200, off200, hz);
  return ret;
}

// return 2 if it decodes to a brand-new message.
// return 1 if it decodes but we've already seen it,
//   perhaps in a different pass.
// return 0 if we could not decode.
int
one_iter(const std::vector<double> &samples200, int best_off, double hz_for_cb)
{
  std::vector<Strength> strengths =
    search_both(samples200,
                25, second_hz_fracs, second_hz_win,
                best_off, second_off_fracs, second_off_win * 32);
  //
  // sort strongest-first.
  //
  std::sort(strengths.begin(), strengths.end(),
            [](const Strength &a, const Strength &b) -> bool
            { return a.strength_ > b.strength_; } );

  for(int i = 0; i < strengths.size() && i < second_count; i++){
    double hz = strengths[i].hz_;
    int off = strengths[i].off_;
    int ret = one_iter1(samples200, off, hz, hz_for_cb, hz_for_cb);
    if(ret > 0){
      return ret;
    }
  }

  return 0;
}

//
// estimate SNR, yielding numbers vaguely similar to WSJT-X.
// m79 is a 79x8 complex FFT output.
//
double
guess_snr(const ffts_t &m79)
{
  int costas[] = { 3, 1, 4, 0, 6, 5, 2 };
  double noises = 0;
  double signals = 0;

  for(int i = 0; i < 7; i++){
    signals += std::abs(m79[i][costas[i]]);
    signals += std::abs(m79[36+i][costas[i]]);
    signals += std::abs(m79[72+i][costas[i]]);
    noises += std::abs(m79[i][(costas[i]+4)%8]);
    noises += std::abs(m79[36+i][(costas[i]+4)%8]);
    noises += std::abs(m79[72+i][(costas[i]+4)%8]);
  }

  for(int i = 0; i < 79; i++){
    if(i < 7 || (i >= 36 && i < 36+7) || (i >= 72 && i < 72+7))
      continue;
    std::vector<double> v(8);
    for(int j = 0; j < 8; j++){
      v[j] = std::abs(m79[i][j]);
    }
    std::sort(v.begin(), v.end());
    signals += v[7]; // strongest tone, probably the signal
    noises += (v[2]+v[3]+v[4])/3;
  }

  noises /= 79;
  signals /= 79;

  noises *= noises; // square yields power
  signals *= signals;

  double raw = signals / noises;
  raw -= 1; // turn (s+n)/n into s/n
  if(raw < 0.1)
    raw = 0.1;
  raw /= (2500.0 / 2.7); // 2.7 hz noise b/w -> 2500 hz b/w
  double snr = 10 * log10(raw);
  snr += 5;
  snr *= 1.4;
  return snr;
}

//
// the signal is at 25 hz in samples200.
// 
// return 2 if it decodes to a brand-new message.
// return 1 if it decodes but we've already seen it,
//   perhaps in a different pass.
// return 0 if we could not decode.
//
int
one_iter1(const std::vector<double> &samples200,
          int best_off, double best_hz,
          double hz0_for_cb, double hz1_for_cb)
{
  // mini 79x8 FFT.
  ffts_t m79 = extract(samples200, best_hz, best_off);

  double ll174[174];
  prepare_soft(m79, ll174, best_off);


  int ret = try_decode(ll174, samples200, best_hz, best_off,
                       hz0_for_cb, hz1_for_cb, 1, "", m79);
  if(ret)
    return ret;

  if(use_hints){
    for(int hi = 0; hi < (int)hints1_.size(); hi++){
      int h = hints1_[hi]; // 28-bit number, goes in ll174 0..28
      if(use_hints == 2 && h != 2){
        // just CQ
        continue;
      }
      double n174[174];
      for(int i = 0; i < 174; i++){
        if(i < 28){
          int bit = h & (1 << 27);
          if(bit){
            n174[i] = -4.97;
          } else {
            n174[i] = 4.97;
          }
          h <<= 1;
        } else {
          n174[i] = ll174[i];
        }
      }
      int ret = try_decode(n174, samples200, best_hz, best_off,
                           hz0_for_cb, hz1_for_cb, 0, "hint1", m79);
      if(ret){
        return ret;
      }
    }
  }
  
  if(use_hints == 1){
    for(int hi = 0; hi < (int)hints2_.size(); hi++){
      int h = hints2_[hi]; // 28-bit number, goes in ll174 29:29+28
      double n174[174];
      for(int i = 0; i < 174; i++){
        if(i >= 29 && i < 29+28){
          int bit = h & (1 << 27);
          if(bit){
            n174[i] = -4.97;
          } else {
            n174[i] = 4.97;
          }
          h <<= 1;
        } else {
          n174[i] = ll174[i];
        }
      }
      int ret = try_decode(n174, samples200, best_hz, best_off,
                           hz0_for_cb, hz1_for_cb, 0, "hint2", m79);
      if(ret){
        return ret;
      }
    }
  }

  return 0;
}

//
// subtract a corrected decoded signal from nsamples_,
// perhaps revealing a weaker signal underneath,
// to be decoded in a subsequent pass.
// re79 is the corrected symbol numbers, as sent over the air.
//
// just zeros out the relevant bin.
// XXX surrounding median amplitude.
// XXX leakage into neighboring bins.
//
void
subtract(const std::vector<int> re79,
         double hz0,
         double hz1,
         double off_sec)
{
  int block = blocksize(rate_);
  double bin_hz = rate_ / (double) block;
  int off0 = off_sec * rate_;

  double mhz = (hz0 + hz1) / 2.0;
  int bin0 = round(mhz / bin_hz);

  // move nsamples so that signal is centered in bin0.
  double diff0 = (bin0 * bin_hz) - hz0;
  double diff1 = (bin0 * bin_hz) - hz1;
  std::vector<double> moved = hilbert_shift(nsamples_, diff0, diff1, rate_);

  ffts_t bins = ffts(moved, off0, block, "subtract1");

  if(bin0 + 8 > bins[0].size())
    return;
  if(bins.size() < 79)
    return;

  // for each symbol, median amplitude of nearby tones.
  std::vector<double> nearby_amp(79);
  if(sub_amp_win > 0){
    for(int si = 0; si < 79; si++){
      std::vector<double> v;
      for(int i = -sub_amp_win; i <= sub_amp_win; i++){
        if(si+i >= 0 && si+i < 79){
          double x = std::abs(bins[si+i][bin0+re79[si+i]]);
          v.push_back(x);
        }
      }
#if 1
      std::sort(v.begin(), v.end());
      nearby_amp[si] = v[v.size() / 2];
#else
      double sum = 0;
      for(int i = 0; i < v.size(); i++)
        sum += v[i];
      nearby_amp[si] = sum / v.size();
#endif
    }
  }

  // for each symbol, median phase of nearby symbols.
  std::vector<double> nearby_phase(79);
  if(fancy_subtract && sub_phase_win > 0){
    std::vector<double> ph(79);
    for(int si = 0; si < 79; si++){
      ph[si] = std::arg(bins[si][bin0+re79[si]]); // -pi .. pi
    }
    
    for(int si = 0; si < 79; si++){
      std::vector<double> v;
      int nhi = 0;
      int nlo = 0;
      for(int i = -sub_phase_win; i <= sub_phase_win; i++){
        if(si+i >= 0 && si+i < 79){
          double x = ph[si+i];
          v.push_back(x);
          if(x < 0)
            nlo++;
          else
            nhi++;
        }
      }

      // try to fix -pi / pi wrap-around.
      int n = v.size();
      if(nhi > nlo){
        for(int i = 0; i < n; i++){
          if(v[i] < -sub_wrap_thresh){
            v[i] += 2*M_PI;
          }
        }
      } else if(nlo > nhi) {
        for(int i = 0; i < n; i++){
          if(v[i] > sub_wrap_thresh){
            v[i] -= 2*M_PI;
          }
        }
      }
      
#if 1
      std::sort(v.begin(), v.end());
      nearby_phase[si] = v[n/2];
#else
      double sum = 0;
      double w = 0;
      for(int i = 0; i < v.size(); i++){
        sum += v[i];
        w += 1;
      }
      nearby_phase[si] = sum / w;
#endif
    }
  }

  for(int si = 0; si < 79; si++){
    int sym = bin0 + re79[si];

    if(fancy_subtract){
      std::complex<double> c = bins[si][sym];

      double phase;
      if(sub_phase_win > 0){
        phase = nearby_phase[si];
      } else {
        phase = std::arg(c);
      }

      double amp;
      if(sub_amp_win > 0){
        amp = nearby_amp[si];
      } else {
        amp = std::abs(c);
      }

      // amp = std::min(amp, std::abs(c));

      // we know what the symbol was, and we have guesses for
      // the signal's phase and amplitude, so subtract.
      
      // subtract in the time domain, since it's pretty involved
      // to do it correctly in the frequency domain. the point is
      // to subtract the likely signal, but to leave behind
      // stuff at that tone that's not the same phase as the signal.
      
      // FFT multiplies magnitudes by number of bins,
      // or half the number of samples.
      amp /= (block / 2);
      
      double theta = phase;
      
      double tonehz = 6.25 * (bin0 + re79[si]);
      
      for(int jj = 0; jj < block; jj++){
        double x = amp * cos(theta);
        int iii = off0 + block*si + jj;
        moved[iii] -= x;
        theta += 2 * M_PI / (rate_ / tonehz);
      }
    } else {
      if(sub_amp_win > 0){
        double aa = std::abs(bins[si][sym]);
        
        double amp = nearby_amp[si];
        if(amp > aa)
          amp = aa;
        
        if(aa > 0.0){
          bins[si][sym] /= aa;
          bins[si][sym] *= (aa - amp);
        }
      } else {
        bins[si][sym] = 0;
      }

      // to get the original amplitude scale, you have to divide
      // by blocksize before calling ifft.
      std::vector<std::complex<double>> bb(bins[si].size());
      for(int i = 0; i < bb.size(); i++)
        bb[i] = bins[si][i] / (double)block;
      
      std::vector<double> ss = one_ifft(bb, "subtract2");
      assert(ss.size() == block);
      for(int jj = 0; jj < block; jj++){
        moved[off0 + block*si + jj] = ss[jj];
      }
    }
  }

  nsamples_ = hilbert_shift(moved, -diff0, -diff1, rate_);
}

//
// decode, give to callback, and subtract.
//
// return 2 if it decodes to a brand-new message.
// return 1 if it decodes but we've already seen it,
//   perhaps in a different pass.
// return 0 if we could not decode.
//
int
try_decode(double ll174[174], const std::vector<double> &samples200,
           double best_hz, int best_off, double hz0_for_cb, double hz1_for_cb,
           int use_osd, const char *comment1,
           const ffts_t &m79)
{
  int a174[174];
  std::string comment(comment1);

  if(decode(ll174, a174, use_osd, comment)){
    // a174 is 91 bits of plain message, 83 bits of parity.

    // reconstruct correct 79 symbols from LDPC output.
    std::vector<int> re79 = recode(a174);

    // fine-tune offset and hz for better subtraction.
    double best_strength = 0;
    search_both_known(samples200, re79,
                      best_hz, third_hz_inc, third_hz_win,
                      best_off, third_off_inc, third_off_win,
                      best_hz, best_off, best_strength);

    double off_sec = best_off / 200.0;

    // hz0_for_cb and corrected_hz* refers to samples_,
    // so that's what we want for subtraction.
    // but we down-shifted what Python gave us by down_hz_,
    // so we also need to add that for the callback.
    
    double corrected_hz0 = hz0_for_cb + (best_hz - 25.0);
    double corrected_hz1 = hz1_for_cb + (best_hz - 25.0);

    double snr = guess_snr(m79);
    
    if(cb_ != 0){
      cb_mu_.lock();
      int ret = cb_(a174, corrected_hz0 + down_hz_, corrected_hz1 + down_hz_,
                    off_sec, comment.c_str(), snr, pass_);
      cb_mu_.unlock();
      if(ret == 2){
        // a new decode.
        // subtract from nsamples_.
        subtract(re79, corrected_hz0, corrected_hz1, off_sec);
      }
      return ret;
    }
    return 1;
  } else {
    return 0;
  }
}

//
// given 174 bits corrected by LDPC, work
// backwards to the symbols that must have
// been sent.
// used to help ensure that subtraction subtracts
// at the right place.
//
std::vector<int>
recode(int a174[])
{
  int i174 = 0;
  int costas[] = { 3, 1, 4, 0, 6, 5, 2 };
  std::vector<int> out79;
  for(int i79 = 0; i79 < 79; i79++){
    if(i79 < 7){
      out79.push_back(costas[i79]);
    } else if(i79 >= 36 && i79 < 36+7){
      out79.push_back(costas[i79-36]);
    } else if(i79 >= 72){
      out79.push_back(costas[i79-72]);
    } else {
      int sym = (a174[i174+0] << 2) | (a174[i174+1] << 1) | (a174[i174+2] << 0);
      i174 += 3;
      // gray code
      int map[] = { 0, 1, 3, 2, 5, 6, 4, 7 };
      sym = map[sym];
      out79.push_back(sym);
    }
  }
  assert(out79.size() == 79);
  assert(i174 == 174);
  return out79;
}

};

std::mutex FT8::cb_mu_;

//
// Python calls these.
//
void
entry(double xsamples[], int nsamples, int start, int rate,
      double min_hz,
      double max_hz,
      int hints1[],
      int hints2[],
      double time_left, double total_time_left, cb_t cb)
{
  double t0 = now();
  double deadline = t0 + time_left;
  double final_deadline = t0 + total_time_left;

  std::vector<double> samples(nsamples);
  for(int i = 0; i < nsamples; i++){
    samples[i] = xsamples[i];
  }

  if(min_hz < 0){
    min_hz = 0;
  }
  if(max_hz > rate/2){
    max_hz = rate/2;
  }
  double per = (max_hz - min_hz) / nthreads;

  std::vector<FT8 *> thv;

  for(int i = 0; i < nthreads; i++){
    double hz0 = min_hz + i * per;
    if(i > 0 || overlap_edges)
      hz0 -= overlap;
    
    double hz1 = min_hz + (i + 1) * per;
    if(i != nthreads-1 || overlap_edges)
      hz1 += overlap;

    FT8 *ft8 = new FT8(samples,
                       hz0, hz1,
                       start, rate,
                       hints1, hints2,
                       deadline, final_deadline, cb);

    ft8->th_ = new std::thread( [ ft8 ] () { ft8->go(); } );
    thv.push_back(ft8);
  }

  for(int i = 0; i < thv.size(); i++){
    thv[i]->th_->join();
    delete thv[i]->th_;
    delete thv[i];
  }
}

double
set(char *param, char *val)
{
  struct sss {
    const char *name;
    void *addr;
    int type; // 0 int, 1 double
  };
  struct sss params[] =
    {
     { "snr_win", &snr_win, 0 },
     { "snr_how", &snr_how, 0 },
     { "soft_ranges", &soft_ranges, 0 },
     { "best_in_noise", &best_in_noise, 0 },
     { "ldpc_iters", &ldpc_iters, 0 },
     { "shoulder200", &shoulder200, 1 },
     { "shoulder200_extra", &shoulder200_extra, 1 },
     { "second_hz_fracs", &second_hz_fracs, 0 },
     { "second_hz_win", &second_hz_win, 1 },
     { "second_off_fracs", &second_off_fracs, 0 },
     { "second_off_win", &second_off_win, 1 },
     { "third_hz_inc", &third_hz_inc, 1 },
     { "third_hz_win", &third_hz_win, 1 },
     { "third_off_inc", &third_off_inc, 0 },
     { "third_off_win", &third_off_win, 0 },
     { "log_tail", &log_tail, 1 },
     { "log_rate", &log_rate, 1 },
     { "problt_how", &problt_how, 0 },
     { "use_apriori", &use_apriori, 0 },
     { "use_hints", &use_hints, 0 },
     { "win_type", &win_type, 0 },
     { "osd_depth", &osd_depth, 0 },
     { "ncoarse", &ncoarse, 0 },
     { "ncoarse_blocks", &ncoarse_blocks, 0 },
     { "tminus", &tminus, 1 },
     { "tplus", &tplus, 1 },
     { "coarse_off_fracs", &coarse_off_fracs, 0 },
     { "coarse_hz_fracs", &coarse_hz_fracs, 0 },
     { "already_hz", &already_hz, 1 },
     { "nthreads", &nthreads, 0 },
     { "npasses", &npasses, 0 },
     { "overlap", &overlap, 1 },
     { "sub_amp_win", &sub_amp_win, 0 },
     { "sub_phase_win", &sub_phase_win, 0 },
     { "nyquist", &nyquist, 1 },
     { "oddrate", &oddrate, 0 },
     { "osd_ldpc_thresh", &osd_ldpc_thresh, 0 },
     { "pass0_frac", &pass0_frac, 1 },
     { "go_extra", &go_extra, 1 },
     { "reduce_len", &reduce_len, 0 },
     { "reduce_how", &reduce_how, 0 },
     { "reduce_margin", &reduce_margin, 1 },
     { "do_reduce", &do_reduce, 0 },
     { "pass_threshold", &pass_threshold, 0 },
     { "strength_how", &strength_how, 0 },
     { "reduce_shoulder", &reduce_shoulder, 1 },
     { "overlap_edges", &overlap_edges, 0 },
     { "coarse_strength_how", &coarse_strength_how, 0 },
     { "coarse_all", &coarse_all, 1 },
     { "second_count", &second_count, 0 },
     { "soft_phase_win", &soft_phase_win, 0 },
     { "sub_wrap_thresh", &sub_wrap_thresh, 1 },
     { "phase_wrap_thresh", &phase_wrap_thresh, 1 },
     { "fancy_subtract", &fancy_subtract, 0 },
    };
  int nparams = sizeof(params) / sizeof(params[0]);

  for(int i = 0; i < nparams; i++){
    if(strcmp(param, params[i].name) == 0){
      if(val[0]){
        if(params[i].type == 0){
          *(int*)params[i].addr = round(atof(val));
        } else if(params[i].type == 1){
          *(double*)params[i].addr = atof(val);
        }
      }
      if(params[i].type == 0){
        return *(int*)params[i].addr;
      } else if(params[i].type == 1){
        return *(double*)params[i].addr;
      } else {
        fprintf(stderr, "weird type %d\n", params[i].type);
        return 0;
      }
    }
  }
  fprintf(stderr, "ft8.cc set(%s, %s) unknown parameter\n", param, val);
  exit(1);
  return 0;
}
