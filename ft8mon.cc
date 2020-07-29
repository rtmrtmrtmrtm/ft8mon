//
// decode FT8 from a sound card
//
// Robert Morris, AB1HL
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <vector>
#include <fftw3.h>
#include <time.h>
#include <string.h>
#include <mutex>
#include <map>
#include <string>
#include "snd.h"
#include "util.h"
#include "unpack.h"

extern "C" {
  typedef int (*cb_t)(int *a91, double hz0, double hz1, double off,
                      const char *, double snr);
  void entry(double xsamples[], int nsamples, int start, int rate,
             double min_hz,
             double max_hz,
             int hints1[], int hints2[], double time_left,
             double total_time_left, cb_t cb);
};

std::mutex cycle_mu;
volatile int cycle_count;
time_t saved_cycle_start;
std::map<std::string,bool> cycle_already;

//
// a91 is 91 bits -- 77 plus the 14-bit CRC.
//
int
hcb(int *a91, double hz0, double hz1, double off,
    const char *comment, double snr)
{
  std::string msg = unpack(a91);

  cycle_mu.lock();

  if(cycle_already.count(msg) > 0){
    // already decoded this message on this cycle
    cycle_mu.unlock();
    return 1; // 1 => already seen, don't subtract.
  }

  cycle_already[msg] = true;
  cycle_count += 1;

  cycle_mu.unlock();

  struct tm result;
  gmtime_r(&saved_cycle_start, &result);

  printf("%02d%02d%02d %3d %4.1f %4d %s\n",
         result.tm_hour,
         result.tm_min,
         result.tm_sec,
         (int)snr,
         off - 1.5,
         (int)hz0,
         msg.c_str());
  fflush(stdout);
  
  return 2; // 2 => new decode, do subtract.
}

void
usage()
{
  fprintf(stderr, "Usage: ft8mon -card card channel\n");
#ifdef AIRSPYHF
  fprintf(stderr, "       ft8mon -card airspy mhz\n");
#endif
  fprintf(stderr, "       ft8mon -levels card channel\n");
  fprintf(stderr, "       ft8mon -list\n");
  fprintf(stderr, "       ft8mon -file xxx.wav\n");
  exit(1);
}

int
main(int argc, char *argv[])
{
  int hints[2] = { 2, 0 }; // CQ
  double budget = 5; // compute for this many seconds per cycle

  // avoid snd card input overruns.
  // it's not clear why there's a problem.
  extern int fftw_type;
  fftw_type = FFTW_ESTIMATE; // rather than FFTW_MEASURE

  extern int nthreads;
  nthreads = 4; // use four cores in parallel
  
  if(argc == 4 && strcmp(argv[1], "-card") == 0){
    SoundIn *sin = SoundIn::open(argv[2], argv[3]);
    sin->start();
    int rate = sin->rate();

    while(1){
      // sleep until 14 seconds into the next 15-second cycle.
      double tt = now();
      long long cycle_start = tt - ((long long)tt % 15);

      if(tt - cycle_start >= 14){
        double ttt_start;
        // get all waiting samples.
        // ttt is UNIX time of samples[0].
        std::vector<double> samples = sin->get(rate * 100, ttt_start);

        double ttt_end = ttt_start + samples.size() * (1.0 / rate);
        cycle_start = ((long long) (ttt_end / 15)) * 15;

        long long nominal = samples.size() - rate * (ttt_end - cycle_start - 0.5);
        if(nominal >= 0){
          printf("decodes: %d\n", cycle_count);

          cycle_mu.lock();
          cycle_count = 0;
          saved_cycle_start = cycle_start; // for hcb() callback
          cycle_already.clear();
          cycle_mu.unlock();

          entry(samples.data(), samples.size(), nominal, rate,
                150,
                2900,
                hints, hints, budget, budget, hcb);
        }

        sleep(2);
      }
      usleep(100 * 1000); // 0.1 seconds
    }
  } else if(argc == 4 && strcmp(argv[1], "-levels") == 0){
    SoundIn *sin = SoundIn::open(argv[2], argv[3]);
    sin->start();
    sin->levels();
  } else if(argc == 3 && strcmp(argv[1], "-file") == 0){
    // the .wav file should start at an even 15-second boundary.
    int rate;
    std::vector<double> s = readwav(argv[2], rate);
    entry(s.data(), s.size(), 0.5 * rate, rate,
          150,
          2900,
          hints, hints, budget, budget, hcb);
  } else if(argc == 2 && strcmp(argv[1], "-list") == 0){
    snd_list();
  } else {
    usage();
  }
}
