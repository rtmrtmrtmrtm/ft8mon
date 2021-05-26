#ifndef snd_h
#define snd_h 1

#include <math.h>
#include <assert.h>
#if defined(USE_AIRSPYHF) || defined(USE_HPSDR) || defined(USE_SDRIP)
#include <liquid/liquid.h>
#endif
#include <complex>
#include <string>
#include <vector>
#include "util.h"
#ifdef USE_AIRSPYHF
#include <airspyhf.h>
#endif
#ifdef USE_HPSDR
#include "hpsdr.h"
#endif
#ifdef USE_SDRIP
#include "sdrip.h"
#endif
#include "cloudsdr.h"

void snd_init();
void snd_list();

class SoundIn {
public:
  virtual void start() = 0;
  virtual int rate() = 0;
  virtual std::vector<double> get(int n, double &t0, int latest) = 0;
  virtual bool has_iq() { return false; } // default
  virtual std::vector<std::complex<double>> get_iq(int n, double &t0, int latest) {
    assert(0);
  }
  void levels();
  virtual int set_freq(int) { return -1; }
  static SoundIn *open(std::string card, std::string chan, int rate);
};

class CardSoundIn : public SoundIn {
 private:
  int card_;
  int chan_;
  int rate_;
  int channels_;
  double dt_; // time difference (seconds) between UNIX and stream time

  // circular buffer
  int n_;
  short *buf_;
  volatile int wi_;
  volatile int ri_;
  double time_; // of most recent sample, buf_[wi_-1]

 public:
  CardSoundIn(int card, int chan, int rate);
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  int rate() { return rate_; }

  static int cb(const void *input,
                void *output,
                unsigned long frameCount,
                const struct PaStreamCallbackTimeInfo *timeInfo,
                unsigned long statusFlags, // PaStreamCallbackFlags
                void *userData);
};

class FileSoundIn : public SoundIn {
private:
  int rate_;
  std::vector<double> samples_;
  int i_;
  double t_; // unix time of samples_[i_];
public:
  FileSoundIn(std::string filename, int rate) {
    samples_ = readwav(filename.c_str(), rate_);
    if(rate != -1 && rate != rate_){
      fprintf(stderr, "FileSoundIn(%s, %d) but rate %d\n",
              filename.c_str(), rate, rate_);
      exit(1);
    }
    i_ = 0;
    t_ = now();
  }
  int rate() { return rate_; }
  void start() { };
  std::vector<double> get(int n, double &t0, int latest) {
    t0 = t_;

    std::vector<double> v;
    for(int j = 0; j < n; j++){
      if (i_ < (int) samples_.size()){
        v.push_back(samples_[i_]);
        i_ += 1;
        t_ += 1 / (double) rate_;
      } else {
        break;
      }
    }

    if(v.size() == 0){
      exit(0); // XXX
    }
        
    return v;
  }
};

#ifdef USE_AIRSPYHF
class AirspySoundIn : public SoundIn {
 private:
  struct airspyhf_device* device_;
  unsigned long long serial_;
  unsigned int hz_;
  int air_rate_; // 192000
  int rate_; // 12000
  double time_; // of most recent sample, buf_[wi_-1]
  long long count_; // of samples, for rate reduction.
  char hostname_[64];

  firfilt_crcf filter_;

  // circular buffer
  int n_;
  std::complex<double> *buf_;
  volatile int wi_;
  volatile int ri_;

  unsigned long long get_serial();

 public:
  AirspySoundIn(std::string chan, int rate);
  ~AirspySoundIn() { firfilt_crcf_destroy(filter_); }
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  int rate() { return rate_; }
  int set_freq(int);

  static int cb1(airspyhf_transfer_t *);
  int cb2(airspyhf_transfer_t *);
};
#endif

#ifdef USE_HPSDR
class HPSDRSoundIn : public SoundIn {
private:
  HPSDR *sdr_;
  int unit_;
  int rate_;
  int hz_;

  // buffered input.
  std::vector<std::complex<double>> buf_;
  double time_; // of most recent sample, buf_[wi_-1]

  void absorb();

public:
  std::vector<double> usb(std::vector<std::complex<double>> &v1);

public:
  HPSDRSoundIn(std::string chan, int rate);
  ~HPSDRSoundIn();
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  bool has_iq() { return true; }
  std::vector<std::complex<double>> get_iq(int n, double &t0, int latest);
  int rate() { return rate_; }
  int set_freq(int);
  HPSDR *sdr() { return sdr_; }
};
#endif

class CloudSoundIn : public SoundIn {
 private:
  std::string chan_;
  CloudSDR *sdr_;

 public:
  CloudSoundIn(std::string chan, int rate);
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  int rate() { return 8000; }
  int set_freq(int);
};

#ifdef USE_SDRIP
class SDRIPSoundIn : public SoundIn {
 private:
  std::string chan_;
  SDRIP *sdr_;
  int rate_; // desired rate
  std::vector<std::complex<double>> buf_;
  double time_; // of most recent sample, buf_[wi_-1]

  std::vector<double> usb(std::vector<std::complex<double>> &v1);
  void absorb();

 public:
  SDRIPSoundIn(std::string chan, int rate);
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  bool has_iq() { return true; }
  int rate() { return rate_; }
  int set_freq(int);
  std::vector<std::complex<double>> get_iq(int n, double &t0, int latest);
};
#endif

class SoundOut {
public:
  virtual void start() = 0;
  virtual int rate() = 0;
  virtual void write(const std::vector<double> &v) = 0;
  virtual void set_freq(int hz) { }
  static SoundOut *open(const std::string card, const std::string chan, int rate);
};

class CardSoundOut : public SoundOut{
 private:
  int card_;
  int rate_;
  void *str_; // PaStream

 public:
  CardSoundOut(int card, int rate);
  void start();
  int rate() { return rate_; }
  void write(const std::vector<short int> &);
  void write(const std::vector<double> &);
};

#ifdef USE_HPSDR
class HPSDRSoundOut : public SoundOut {
private:
  int rate_;
  HPSDR *sdr_;
  int hz_; // carrier e.g. 14070000
  
public:
  HPSDRSoundOut(const std::string chan, int rate);
  void start();
  int rate() { return rate_; }
  void write(const std::vector<double> &);
  void set_freq(int hz) { hz_ = hz; }
};
#endif

#endif
