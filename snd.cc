#include "snd.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>
#include "util.h"
#include "fft.h"
#include <portaudio.h>
#include <ctype.h>

void
snd_init()
{
  static int inited = 0;

  if(!inited){
    inited = 1;
    int err = Pa_Initialize();
    if(err != paNoError){
      fprintf(stderr, "Pa_Initialize() failed\n");
      exit(1);
    }
  }
}

//
// print a list of sound devices
//
void
snd_list()
{
  snd_init();
  int n = Pa_GetDeviceCount();
  printf("%d sound devices:\n", n);
  for(int di = 0; di < n; di++){
    const PaDeviceInfo *info = Pa_GetDeviceInfo(di);
    if(info == 0)
      continue;
    printf("%d: %s %d/%d ",
           di,
           info->name,
           info->maxInputChannels,
           info->maxOutputChannels);

    PaStreamParameters ip;
    memset(&ip, 0, sizeof(ip));
    ip.device = di;
    ip.channelCount = 1;
    ip.sampleFormat = paInt16;
    ip.suggestedLatency = 0;
    ip.hostApiSpecificStreamInfo = 0;
    int rates[] = { 6000, 8000, 11025, 12000, 16000, 22050, 44100, 48000, 0 };
    for(int ri = 0; rates[ri]; ri++){
      PaError err = Pa_IsFormatSupported(&ip, 0, rates[ri]);
      if(err == paNoError){
        printf("%d ", rates[ri]);
      }
    }

    // printf("-- %f %f ", info->defaultLowOutputLatency, info->defaultHighOutputLatency);

    printf("\n");
  }

#ifdef AIRSPYHF
  int ndev = airspyhf_list_devices(0, 0);
  if(ndev > 0){
    unsigned long long serials[20];
    if(ndev > 20)
      ndev = 20;
    airspyhf_list_devices(serials, ndev);
    for(int i = 0; i < ndev; i++){
      airspyhf_device_t *dev = 0;
      if(airspyhf_open_sn(&dev, serials[i]) == AIRSPYHF_SUCCESS){
        airspyhf_read_partid_serialno_t read_partid_serialno;
        airspyhf_board_partid_serialno_read(dev, &read_partid_serialno);
        printf("Airspy HF+ S/N %08X%08X\n",
               read_partid_serialno.serial_no[0],
               read_partid_serialno.serial_no[1]);
      } else {
        fprintf(stderr, "could not open airspyhf serial %llu\n", serials[i]);
      }
    }
  }
#endif
}

//
// print avg and peak each second.
//
void
SoundIn::levels()
{
  double max = 0;
  double sum = 0;
  int n = 0;
  double last_t = now();

  while(1){
    double dummy;
    std::vector<double> buf = get(rate(), dummy);
    if(buf.size() == 0)
      usleep(100*1000);
    for(int i = 0; i < buf.size(); i++){
      sum += fabs(buf[i]);
      n += 1;
      if(fabs(buf[i]) > max){
        max = fabs(buf[i]);
      }
      if(n >= rate()){
        printf("avg=%.3f peak=%.3f rate=%.1f\n", sum / n, max, n / (now() - last_t));
        n = 0;
        sum = 0;
        max = 0;
        last_t = now();
      }
    }
  }
}

//
// generic open
//
SoundIn *
SoundIn::open(std::string card, std::string chan)
{
  assert(card.size() > 0);
  SoundIn *sin;

  if(isdigit(card[0])){
    sin = new CardSoundIn(atoi(card.c_str()), atoi(chan.c_str()));
  } else if(card == "file"){
    sin = new FileSoundIn(chan);
#ifdef AIRSPYHF
  } else if(card == "airspy"){
    sin = new AirspySoundIn(chan);
#endif
  } else {
    fprintf(stderr, "SoundIn::open(%s, %s): type not recognized\n", card.c_str(), chan.c_str());
    exit(1);
  }

  return sin;
}
  
int
CardSoundIn::cb(const void *input,
            void *output,
            unsigned long frameCount,
            const struct PaStreamCallbackTimeInfo *timeInfo,
            unsigned long statusFlags, // PaStreamCallbackFlags
            void *userData)
{
  CardSoundIn *sin = (CardSoundIn *) userData;
  const short int *buf = (const short int *) input;

  if(statusFlags != 0){
    // 2 is paInputOverflow
    fprintf(stderr, "CardSoundIn::cb statusFlags 0x%x\n", (int)statusFlags);
    exit(1);
  }

  for(int i = 0; i < frameCount; i++){
    if(((sin->wi_ + 1) % sin->n_) != sin->ri_){
      sin->buf_[sin->wi_] = buf[i*sin->channels_ + sin->chan_];
      sin->wi_ = (sin->wi_ + 1) % sin->n_;
    } else {
      fprintf(stderr, "CardSoundIn::cb buf_ overflow\n");
      exit(1);
      break;
    }
  }

  sin->time_ = timeInfo->inputBufferAdcTime + frameCount * (1.0 / sin->rate_);

  return 0;
}

CardSoundIn::CardSoundIn(int card, int chan)
{
  card_ = card;
  chan_ = chan;
  assert(chan_ >= 0 && chan_ <= 1);
  time_ = -1;
}

//
// read a bunch of recent sound samples.
// read up to n samples, no more.
// return immediately with whatever samples exist,
// perhaps fewer than n.
// return UNIX time of first sample in t0.
//
std::vector<double>
CardSoundIn::get(int n, double &t0)
{
  std::vector<double> v;

  if(time_ < 0 && wi_ == ri_){
    // no input has ever arrived.
    t0 = -1;
    return v;
  }

  // calculate time of first sample in buf_.
  // XXX there's a race here with cb().
  t0 = time_ + dt_; // time of last sample in buf_.
  if(wi_ > ri_){
    t0 -= (wi_ - ri_) * (1.0 / rate_);
  } else {
    t0 -= ((wi_ + n_) - ri_) * (1.0 / rate_);
  }

  while(v.size() < n){
    if(ri_ == wi_){
      break;
    }
    short x = buf_[ri_];
    v.push_back(x / 32767.0);
    ri_ = (ri_ + 1) % n_;
  }

  return v;
}

void
CardSoundIn::start()
{
  snd_init();

#ifdef __linux__
  // RIGblaster only supports 44100 and 48000.
  rate_ = 48000;
#else
  rate_ = 6000; // on some systems this may need to be 12000
#endif

#ifdef __FreeBSD__
  // must read both, otherwise FreeBSD mixes them.
  channels_ = 2;
#else
  if(chan_ == 0){
    channels_ = 1;
  } else {
    channels_ = 2;
  }
#endif

  PaStreamParameters ip;
  memset(&ip, 0, sizeof(ip));
  ip.device = card_;
  ip.channelCount = channels_;
  ip.sampleFormat = paInt16;
  // don't set latency to zero; this causes problems on Linux.
  ip.suggestedLatency = Pa_GetDeviceInfo(card_)->defaultLowInputLatency;
  ip.hostApiSpecificStreamInfo = 0;
  
  PaStream *str = 0;
  PaError err = Pa_OpenStream(&str,
                              &ip,
                              0,
                              rate_,
#ifdef __FreeBSD__
                              128, // framesPerBuffer
#else
                              0, // framesPerBuffer
#endif
                              0,
                              cb,
                              (void*) this);
  if(err != paNoError){
    fprintf(stderr, "Pa_OpenStream(card=%d,rate=%d) failed for input: %s\n",
            card_, rate_, Pa_GetErrorText(err));
    exit(1);
  }

  // allocate a 30-second circular buffer
  n_ = rate_ * 30;
  buf_ = (short *) malloc(sizeof(short) * n_);
  assert(buf_);
  wi_ = 0;
  ri_ = 0;

  err = Pa_StartStream(str);
  if(err != paNoError){
    fprintf(stderr, "Pa_StartStream failed\n");
    exit(1);
  }

  dt_ = now() - Pa_GetStreamTime(str);

}

#ifdef AIRSPYHF
//
// chan argument is megahertz.
//
AirspySoundIn::AirspySoundIn(std::string chan)
{
  device_ = 0;
  hz_ = 1000000.0 * atof(chan.c_str());
  air_rate_ = 192 * 1000;;
  rate_ = 12000;
  total_ = 0;
  start_time_ = 0;

  if( airspyhf_open(&device_) != AIRSPYHF_SUCCESS ) {
    fprintf(stderr, "airspyhf_open() failed\n");
    exit(1);
  }

  if (airspyhf_set_samplerate(device_, air_rate_) != AIRSPYHF_SUCCESS) {
    fprintf(stderr, "airspyhf_set_samplerate(%d) failed\n", air_rate_);
    exit(1);
  }

  // allocate a 30-second circular buffer
  n_ = rate_ * 30;
  buf_ = (std::complex<double> *) malloc(sizeof(std::complex<double>) * n_);
  assert(buf_);
  wi_ = 0;
  ri_ = 0;
  time_ = -1;

  // Liquid DSP filter + resampler to convert 192000 to 12000.
  // firdecim?
  // iirdecim?
  // msresamp?
  // msresamp2?
  // resamp?
  // resamp2?
  int h_len = estimate_req_filter_len(0.01, 60.0);
  float h[h_len];
  double cutoff = (rate_ / (double) air_rate_) / 2.0;
  liquid_firdes_kaiser(h_len,
                       cutoff,
                       60.0,
                       0.0,
                       h);
  assert((air_rate_ % rate_) == 0);
  filter_ = firfilt_crcf_create(h, h_len);
}

void
AirspySoundIn::start()
{
  if( airspyhf_start(device_, cb1, this) != AIRSPYHF_SUCCESS ) {
    fprintf(stderr, "airspyhf_start() failed.\n");
    exit(1);
  }

  if( airspyhf_set_freq(device_, hz_) != AIRSPYHF_SUCCESS ) {
    fprintf(stderr, "airspyhf_set_freq(%d) failed.\n", hz_);
    exit(1);
  }
}

int
AirspySoundIn::cb1(airspyhf_transfer_t *transfer)
{
  AirspySoundIn *sin = (AirspySoundIn *) transfer->ctx;
  return sin->cb2(transfer);
}

int
AirspySoundIn::cb2(airspyhf_transfer_t *transfer)
{
  // transfer->sample_count
  // transfer->samples
  // transfer->ctx
  // transfer->dropped_samples
  // each sample is an I/Q pair of 32-bit floats
  if(total_ == 0){
    start_time_ = now() - (1.0 / air_rate_) * transfer->sample_count;
  }

  if(transfer->dropped_samples){
    fprintf(stderr, "dropped_samples %d\n", (int)transfer->dropped_samples);
  }

  // I/Q
  // low-pass-filter
  //   maybe this: https://liquidsdr.org/doc/msresamp/

  airspyhf_complex_float_t *buf = transfer->samples;
  for(int i = 0; i < transfer->sample_count; i++){
    // low-pass filter, preparatory to rate reduction.
    liquid_float_complex x, y;
    x.real = buf[i].re;
    x.imag = buf[i].im;
    firfilt_crcf_push(filter_, x);
    firfilt_crcf_execute(filter_, &y);

    if((total_ % (air_rate_ / rate_)) == 0){
      if(((wi_ + 1) % n_) != ri_){
        // XXX what is the range of buf[i].re? 0..1?
        buf_[wi_] = std::complex<double>(y.real, y.imag);
        wi_ = (wi_ + 1) % n_;
      } else {
        fprintf(stderr, "CardSoundIn::cb buf_ overflow\n");
        exit(1); // XXX
        break;
      }
    }
    total_ += 1;
  }

  time_ = start_time_ + total_ * (1.0 / air_rate_);

  return 0;
}

std::vector<double>
vreal(std::vector<std::complex<double>> a)
{
  std::vector<double> b(a.size());
  for(int i = 0; i < a.size(); i++){
    b[i] = a[i].real();
  }
  return b;
}

std::vector<double>
vimag(std::vector<std::complex<double>> a)
{
  std::vector<double> b(a.size());
  for(int i = 0; i < a.size(); i++){
    b[i] = a[i].imag();
  }
  return b;
}

//
// convert I/Q to USB via the phasing method.
// uses FFTs over the whole of a[], so it's slow.
// and the results are crummy at the start and end.
//
std::vector<double>
iq2usb(std::vector<std::complex<double>> a)
{
  std::vector<double> ii = vreal(analytic(vreal(a)));
  std::vector<double> qq = vimag(analytic(vimag(a)));
  std::vector<double> ssb(ii.size());
  for(int i = 0; i < ii.size(); i++){
    ssb[i] = ii[i] - qq[i];
  }
  return ssb;
}

std::vector<double>
AirspySoundIn::get(int n, double &t0)
{
  std::vector<double> nothing;

  if(time_ < 0 && wi_ == ri_){
    // no input has ever arrived.
    t0 = -1;
    return nothing;
  }

  // calculate time of first sample in buf_.
  // XXX there's a race here with cb().
  t0 = time_; // time of last sample in buf_.
  if(wi_ > ri_){
    t0 -= (wi_ - ri_) * (1.0 / rate_);
  } else {
    t0 -= ((wi_ + n_) - ri_) * (1.0 / rate_);
  }

  std::vector<std::complex<double>> v1;
  while(v1.size() < n){
    if(ri_ == wi_){
      break;
    }
    v1.push_back(buf_[ri_]);
    ri_ = (ri_ + 1) % n_;
  }

  if(v1.size() < 2){
    // analytic() demands more than one sample.
    return vreal(v1);
  } else {
    std::vector<double> v2 = iq2usb(v1);
    return v2;
  }
}
#endif

SoundOut::SoundOut(int card)
{
  card_ = card;
}

void
SoundOut::start()
{
  snd_init();

#ifdef __linux__
  // RIGblaster only supports 44100 and 48000.
  rate_ = 48000;
#else
  rate_ = 8000;
#endif

  PaStreamParameters op;
  memset(&op, 0, sizeof(op));
  op.device = card_;
  op.channelCount = 1;
  op.sampleFormat = paInt16;
  // latency must be the same as for input card.
  op.suggestedLatency = Pa_GetDeviceInfo(card_)->defaultLowInputLatency;
  op.hostApiSpecificStreamInfo = 0;
  
  str_ = 0;
  PaError err = Pa_OpenStream(&str_,
                              0,
                              &op,
                              rate_,
                              0, // framesPerBuffer
                              0,
                              0,
                              (void*) 0);
  if(err != paNoError){
    fprintf(stderr, "Pa_OpenStream(card=%d,rate=%d) failed for output: %s\n",
            card_, rate_, Pa_GetErrorText(err));
    exit(1);
  }

  err = Pa_StartStream(str_);
  if(err != paNoError){
    fprintf(stderr, "Pa_StartStream failed\n");
    exit(1);
  }

}

//
// Pa_WriteStream may block, but it does have some internal buffering.
// the amount seems to be controlled by suggestedLatency.
//
void
SoundOut::write(const std::vector<short int> &v)
{
  PaError err = Pa_WriteStream(str_, v.data(), v.size());
  if(err != paNoError && err != paOutputUnderflowed){
    fprintf(stderr, "Pa_WriteStream failed %d %s\n", err, Pa_GetErrorText(err));
    exit(1);
  }
}

void
SoundOut::write(const std::vector<double> &v)
{
  std::vector<short int> vv(v.size());
  for(int i = 0; i < v.size(); i++){
    if(v[i] > 1.0){
      fprintf(stderr, "SoundOut::write() oops %f\n", v[i]);
    }
    vv[i] = v[i] * 16380;
  }
  write(vv);
}
