//
// main() routine to decode FT8 from a sound card.
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

time_t saved_cycle_start;

//
// a91 is 91 bits -- 77 plus the 14-bit CRC.
//
int
hcb(int *a91, double hz0, double hz1, double off,
    const char *comment, double snr)
{
  std::string msg = unpack(a91);

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
  
  return 2; // indicate it's a good-looking new decode.
}

void
usage()
{
  fprintf(stderr, "Usage: ft8mon -card card-num channel\n");
  snd_list();
  exit(1);
}

int
main(int argc, char *argv[])
{
  // avoid snd card input overruns.
  // it's not clear why there's a problem.
  extern int fftw_type;
  fftw_type = FFTW_ESTIMATE; // rather than FFTW_MEASURE

  extern int nthreads;
  nthreads = 4;
  
  if(argc == 4 && strcmp(argv[1], "-card") == 0){
    CardSoundIn *sin = new CardSoundIn(atoi(argv[2]), atoi(argv[3]));
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
          int hints[2] = { 2, 0 };
          saved_cycle_start = cycle_start; // for hcb() callback
          entry(samples.data(), samples.size(), nominal, rate,
                150,
                2900,
                hints, hints, 4, 6, hcb);
        }

        sleep(2);
      }
      usleep(100 * 1000); // 0.1 seconds
    }
  } else {
    usage();
  }
}
