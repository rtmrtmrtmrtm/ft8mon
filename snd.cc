#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>
#include "snd.h"
#include "util.h"
#include <portaudio.h>

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
}

//
// print avg and peak each second.
//
void
levels(int card, int chan)
{
  snd_init();
  CardSoundIn *sin = new CardSoundIn(card, chan);
  sin->start();

  double max = 0;
  double sum = 0;
  int n = 0;

  while(1){
    double dummy;
    std::vector<double> buf = sin->get(1024, dummy);
    if(buf.size() == 0)
      usleep(100*1000);
    for(int i = 0; i < buf.size(); i++){
      sum += fabs(buf[i]);
      n += 1;
      if(fabs(buf[i]) > max){
        max = fabs(buf[i]);
      }
      if(n >= sin->rate()){
        printf("%f %f\n", sum / n, max);
        n = 0;
        sum = 0;
        max = 0;
      }
    }
  }
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
    exit(1); // XXX
  }

  for(int i = 0; i < frameCount; i++){
    if(((sin->wi_ + 1) % sin->n_) != sin->ri_){
      sin->buf_[sin->wi_] = buf[i*sin->channels_ + sin->chan_];
      sin->wi_ = (sin->wi_ + 1) % sin->n_;
    } else {
      fprintf(stderr, "CardSoundIn::cb buf_ overflow\n");
      exit(1); // XXX
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
  ip.suggestedLatency = 0;
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
  op.suggestedLatency = 0.1; // avoid output glitches due to underflow.
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
