#ifndef snd_h
#define snd_h 1

#include <math.h>
#include <assert.h>
#ifdef AIRSPYHF
#include <liquid/liquid.h>
#endif
#include <complex>
#include <string>
#include <vector>
#include "util.h"
#ifdef AIRSPYHF
#include <airspyhf.h>
#endif

void snd_init();
void snd_list();

class SoundIn {
public:
  virtual void start() = 0;
  virtual int rate() = 0;
  virtual std::vector<double> get(int n, double &t0, int latest) = 0;
  void levels();
  virtual int set_freq(int) { return -1; }
  static SoundIn *open(std::string card, std::string chan);
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
  CardSoundIn(int card, int chan);
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
  FileSoundIn(std::string filename) {
    samples_ = readwav(filename.c_str(), rate_);
    i_ = 0;
    t_ = now();
  }
  int rate() { return rate_; }
  void start() { };
  std::vector<double> get(int n, double &t0, int latest) {
    t0 = t_;

    std::vector<double> v;
    for(int j = 0; j < n; j++){
      if (i_ < samples_.size()){
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

#ifdef AIRSPYHF
class AirspySoundIn : public SoundIn {
 private:
  struct airspyhf_device* device_;
  unsigned int hz_;
  int air_rate_; // 192000
  int rate_; // 12000
  double time_; // of most recent sample, buf_[wi_-1]
  long long count_; // of samples, for rate reduction.

  firfilt_crcf filter_;

  // circular buffer
  int n_;
  std::complex<double> *buf_;
  volatile int wi_;
  volatile int ri_;

 public:
  AirspySoundIn(std::string chan);
  ~AirspySoundIn() { firfilt_crcf_destroy(filter_); }
  void start();
  std::vector<double> get(int n, double &t0, int latest);
  int rate() { return rate_; }
  int set_freq(int);

  static int cb1(airspyhf_transfer_t *);
  int cb2(airspyhf_transfer_t *);
};
#endif

class SoundOut {
 private:
  int card_;
  int rate_;
  void *str_; // PaStream

 public:
  SoundOut(int card);
  void start();
  int rate() { return rate_; }
  void write(const std::vector<short int> &);
  void write(const std::vector<double> &);
};

#endif
