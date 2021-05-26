//
// Talk over the ethernet to an OpenHPSDR radio
// using Protocol 2.
//
// Tested with a blue Apache Labs ANAN-7000dle.
//
// Robert Morris, AB1HL.
//

#include "hpsdr.h"
#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <sys/select.h>
#include <vector>
#include <thread>
#include <sndfile.h>
#include <ifaddrs.h>
#include <net/if.h>
#include "fft.h"
#include "util.h"

//
// protocol 2 packet formats from:
// openHPSDR Ethernet Protocol V3.8,
// March 2019, Phil Harman VK6PH
//

//
// the sequence:
// host sends discovery packet to broadcast:1024.
//   radio replies, revealing its IP address.
// host sends setup packets (but radio doesn't reply):
//   general registers -> port 1024
//   ddc registers -> port 1025
//   duc registers -> port 1026
//   high priority -> port 1027 (with "run" bit set)
// radio sends:
//   high priority status from 1025
//   mic in from 1026
//   wideband from 1027
//   DDC0 I&Q rom 1035
//
// host must send a Command & Control packet (C&C) every 100 ms.
//

#define ANGELIA_BOARD 3
#define ORION2_BOARD 5

// append a 32-bit value.
// most significant byte first.
static void
push32(std::vector<char> &v, int x)
{
  v.push_back((x >> 24) & 0xff);
  v.push_back((x >> 16) & 0xff);
  v.push_back((x >> 8) & 0xff);
  v.push_back((x >> 0) & 0xff);
}

// x may be negative.
static void
push24(std::vector<char> &v, int x)
{
  if(x > 0)
    assert(x < (1<<23));
  if(x < 0)
    assert(x >= -(1<<23));
  v.push_back((x >> 16) & 0xff);
  v.push_back((x >> 8) & 0xff);
  v.push_back((x >> 0) & 0xff);
}

static void
push16(std::vector<char> &v, int x)
{
  v.push_back((x >> 8) & 0xff);
  v.push_back((x >> 0) & 0xff);
}

static void
push8(std::vector<char> &v, int x)
{
  v.push_back((x >> 0) & 0xff);
}

static unsigned int
r32(const std::vector<char> &v, int off)
{
  unsigned int x = 0;
  x |= (v[off+0] & 0xff) << 24;
  x |= (v[off+1] & 0xff) << 16;
  x |= (v[off+2] & 0xff) << 8;
  x |= (v[off+3] & 0xff) << 0;
  return x;
}

//
// signed! for I/Q samples.
//
static int
r24(const std::vector<char> &v, int off)
{
  unsigned int x = 0;
  x |= (v[off+0] & 0xff) << 16;
  x |= (v[off+1] & 0xff) << 8;
  x |= (v[off+2] & 0xff) << 0;
  if(x & (1 << 23))
    x |= 0xff000000;
  return x;
}

static unsigned int
r16(const std::vector<char> &v, int off)
{
  unsigned int x = 0;
  x |= (v[off+0] & 0xff) << 8;
  x |= (v[off+1] & 0xff) << 0;
  return x;
}

HPSDR *
HPSDR::open()
{
  static HPSDR *sdr = 0;

  if(sdr == 0){
    sdr = new HPSDR;
    std::thread th( [ ] () { sdr->loop(); } );
    th.detach();
  }

  return sdr;
}

HPSDR::HPSDR()
{
  state_ = None;
  s_ = -1;
  last_discover_ = 0;
  last_general_ = 0;
  last_high_ = 0;
  last_swr_print_ = 0;
  run_ = 0;
  board_type_ = -1;
  do_dither_ = 1;
  do_random_ = 1;
  do_lpf_ = 1;
  do_hpf_ = 1;
  do_attenuate_ = 0;
  do_ext1_ = 0;
  do_both_ = 0;
  step_[0] = 0;
  step_[1] = 0;
  last_clipped_[0] = 0;
  last_clipped_[1] = 0;
  last_unclipped_[0] = 0;
  last_unclipped_[1] = 0;

  tx_active_ = 0;
  tx_seq_ = 0;
  tx_ptt_ = 0;
  tx_power_ = 1;
  tx_hz_ = 0;
  tx_pa_ = 0; // enable power amplifier?
  tx_drive_ = 0;
  tx_drive_override_ = -1;

  for(int i = 0; i < NDDC; i++){
    ddc_[i].active_ = 0;
    ddc_[i].hz_ = 0;
    ddc_[i].sdr_rate_ = 0;
    ddc_[i].rate_ = 0;
    ddc_[i].count_ = 0;
    ddc_[i].expect_seq_ = 0;
  }

  check_config_files();
}

//
// allocate an unused DDC.
//
int
HPSDR::allocate_unit(int rate)
{
  for(int unit = 0; unit < NDDC; unit++){
    if(ddc_[unit].hz_ == 0){
      activate(unit, rate);
      ddc_[unit].active_ = 1;
      configuration_changed_ = 1;
      return unit;
    }
  }
  return -1;
}

void
HPSDR::set_freq(int unit, int hz)
{
  assert(unit >= 0 && unit < NDDC && ddc_[unit].active_);
  ddc_[unit].hz_ = hz;
}

void
HPSDR::activate(int unit, int rate)
{
  assert(unit >= 0 && unit < NDDC);
  
  ddc_[unit].sdr_rate_ = 48000;
  if(rate > 0)
    ddc_[unit].rate_ = rate;
  else
    ddc_[unit].rate_ = ddc_[unit].sdr_rate_;
  assert(ddc_[unit].rate_ <= ddc_[unit].sdr_rate_);
  assert((ddc_[unit].sdr_rate_ % ddc_[unit].rate_) == 0);
  ddc_[unit].count_ = 0;

  if(ddc_[unit].rate_ < ddc_[unit].sdr_rate_){
    // Liquid DSP filter + resampler to convert sdr_rate_ to rate_.
    int h_len = estimate_req_filter_len(0.01, 60.0);
    float h[h_len];
    double cutoff = (ddc_[unit].rate_ / (double) ddc_[unit].sdr_rate_) / 2.0;
    
#if 1
    float bands[4] = {
      0.0, (float)(cutoff * 0.9), // pass-band
      (float)cutoff, 0.5 };     // stop-band
    float des[2] = { 1.0, 0.0 }; // desired response
    float weights[2] = { 1.0, 1.0 };
    liquid_firdespm_wtype wtype[2] = {
      LIQUID_FIRDESPM_EXPWEIGHT,
      LIQUID_FIRDESPM_FLATWEIGHT,
    };
    liquid_firdespm_btype btype = LIQUID_FIRDESPM_BANDPASS;
    firdespm_run(h_len, 2, bands, des, weights, wtype, btype, h);
#else
    cutoff *= 0.9;
    liquid_firdes_kaiser(h_len, cutoff, 60.0, 0.0, h);
#endif
    
    ddc_[unit].filter_ = firfilt_crcf_create(h, h_len);
  }
}

HPSDR::~HPSDR()
{
  if(s_ >= 0)
    close(s_);
}

void
HPSDR::list()
{
  int s = open_socket();
  discover(s); // send a discovery broadcast packet.
  for(int i = 0; i < 20 && wait_socket(s, 100); i++){
    struct sockaddr_in from;
    socklen_t fromlen = sizeof(from);
    char buf[2048];
    int cc = recvfrom(s, buf, sizeof(buf), 0, (sockaddr *) &from, &fromlen);
    if(cc < 0){
      perror("hpsdr: recvfrom");
      exit(1);
    }

    if(cc == 60 && (buf[4] == 0x02 || buf[4] == 0x03)){
      // 0x02 idle, 0x03 busy
      const char *types[] = {
        "Atlas", "Hermes", "Hermes2", "Angelia", "Orion",
        "OrionMkII", "HermesLite" };
      int type = buf[11] & 0xff;
      int protocol = buf[12] & 0xff;
      int firmware = buf[13] & 0xff;
      int ddcs = buf[20] & 0xff; // number of DDCs
        
      printf("HPSDR %s %s %s protocol=%d.%d firmware=%d.%d ddcs=%d\n",
             inet_ntoa(from.sin_addr),
             (type >= 0 && type < 7) ? types[type] : "???",
             buf[4] == 0x02 ? "idle" : "busy",
             protocol / 10,
             protocol % 10,
             firmware / 10,
             firmware % 10,
             ddcs);
    }
  }

  close(s);
}

std::vector<char>
HPSDR::make_discovery_packet()
{
  std::vector<char> v;
    
  push32(v, 0); // Seq #
  v.push_back(0x02); // Command - Discovery
  while(v.size() < 60)
    v.push_back(0);
    
  return v;
}

int
HPSDR::open_socket()
{
  int s = socket(AF_INET, SOCK_DGRAM, 0);
  if(s < 0){
    perror("hpsdr: socket");
    exit(1);
  }
    
  int yes = 1;
  if(setsockopt(s, SOL_SOCKET, SO_BROADCAST, &yes, sizeof(yes)) < 0){
    perror("hpsdr: SO_BROADCAST");
  }

  int aaa = 0;
  socklen_t aaalen = sizeof(aaa);
  getsockopt(s, SOL_SOCKET, SO_RCVBUF, &aaa, &aaalen);
  int bbb = 1024 * 1024;
  if(bbb > aaa){
    if(setsockopt(s, SOL_SOCKET, SO_RCVBUF, &bbb, sizeof(bbb)) < 0){
      perror("hpsdr: SO_RCVBUF");
    }
    // int ccc = 0;
    // socklen_t ccclen  = sizeof(ccc);
    // getsockopt(s, SOL_SOCKET, SO_RCVBUF, &ccc, &ccclen);
    // fprintf(stderr, "RCVBUF %d %d %d\n", aaa, bbb, ccc);
  }

  return s;
}

//
// send out a broadcast discover packet on each network interface.
//
void
HPSDR::discover(int s)
{
  struct ifaddrs *ifap;
  if(getifaddrs(&ifap) < 0){
    perror("hpsdr: getifaddrs");
    exit(1);
  }

  for(struct ifaddrs *ifa = ifap; ifa; ifa = ifa->ifa_next){
    if(ifa->ifa_broadaddr == 0)
      continue;
    if((ifa->ifa_flags & IFF_UP) == 0)
      continue;
    if((ifa->ifa_flags & IFF_BROADCAST) == 0)
      continue;
    if(ifa->ifa_broadaddr->sa_family != AF_INET)
      continue;

    sockaddr_in sin;
    memset(&sin, 0, sizeof(sin));
    sin.sin_family = AF_INET;

    // despite what the protocol 2 document says, it seems like
    // you have to send to the broadcast address, not to
    // a specific IP address.

    sin.sin_addr = ((sockaddr_in *)(ifa->ifa_broadaddr))->sin_addr;

    // port 1024 is the discovery port
    sin.sin_port = htons(1024);
    
    std::vector<char> pkt = make_discovery_packet();
    if(sendto(s, pkt.data(), pkt.size(), 0, (sockaddr *) &sin, sizeof(sin)) < 0){
      perror("hpsdr: sendto");
    }
  }

  freeifaddrs(ifap);
}

//
// if file exists, read a number from it and put it in flag.
// otherwise do nothing.
//
void
check_file(const char *file, int &flag)
{
  FILE *fp = fopen(file, "r");
  if(fp == 0)
    return;
  int xxx = 0;
  if(fscanf(fp, "%d", &xxx) == 1){
    if(flag != xxx)
      fprintf(stderr, "%s changed from %d to %d\n", file, flag, xxx);
    flag = xxx;
  }
  fclose(fp);
}

void
check_file(const char *file, volatile int &flag)
{
  int z = flag;
  check_file(file, z);
  flag = z;
}

void
HPSDR::check_config_files()
{
  check_file("hpsdr-dither", do_dither_);
  check_file("hpsdr-random", do_random_);
  check_file("hpsdr-lpf", do_lpf_);
  check_file("hpsdr-hpf", do_hpf_);
  check_file("hpsdr-attenuate", do_attenuate_);
  check_file("hpsdr-ext1", do_ext1_);
  check_file("hpsdr-pa", tx_pa_);
  check_file("hpsdr-power", tx_power_);
  check_file("hpsdr-drive", tx_drive_override_);
  check_file("hpsdr-both", do_both_);
}

//
// run timers.
//
void
HPSDR::tick()
{
  if(state_ == None){
    s_ = open_socket();
    state_ = Discovering;
    run_ = 1;
  }
  if(state_ == Discovering && now() > last_discover_ + 5){
    discover(s_);
    last_discover_ = now();
  }
  if(state_ == Running && (configuration_changed_ || now() > last_general_ + 5)){
    configuration_changed_ = 0;
    
    if(run_ && (int) ddc_[0].input_.size() > 100*ddc_[0].rate_){
      fprintf(stderr, "HPSDR::tick too much input\n");
      ddc_[0].input_.clear();
    }

    check_config_files();

    send_general();
    send_ddc();
    send_duc();
    last_general_ = now();
  }
  if(state_ == Running && now() > last_high_ + 0.2){
    // send "high priority" packets frequently, to keep
    // the watchdog timer from expiring, and to change
    // frequency promptly if need be.
    send_high();
    last_high_ = now();
  }
}

//
// wait up to ms for input on s_,
// return 1 if s_ is readable, 0 otherwise.
//
int
HPSDR::wait_socket(int s, int ms)
{
  fd_set readfds, writefds, exceptfds;
  FD_ZERO(&readfds);
  FD_ZERO(&writefds);
  FD_ZERO(&exceptfds);
  FD_SET(s, &readfds);
  timeval tv;
  tv.tv_sec = ms / 1000;
  tv.tv_usec = (ms % 1000) * 1000;
  int n = select(s + 1, &readfds, &writefds, &exceptfds, &tv);
  if(n > 0 && FD_ISSET(s, &readfds)){
    return 1;
  } else {
    return 0;
  }
}

void
HPSDR::process()
{
  struct sockaddr_in from;
  socklen_t fromlen = sizeof(from);
  char buf[2048];
  int cc = recvfrom(s_, buf, sizeof(buf), 0, (sockaddr *) &from, &fromlen);
  if(cc < 0){
    perror("hpsdr: recvfrom");
    exit(1);
  }
  process_packet(std::vector<char>(buf, buf+cc), from);
}

void
HPSDR::process_discovery_response(const std::vector<char> &buf,
                                  const sockaddr_in &from)
{
  radio_sin_ = from;
  board_type_ = buf[11] & 0xff; // 3 = Angelia, 5 = Orion Mk II
}

//
// send general setup packet.
//
void
HPSDR::send_general()
{
  std::vector<char> v;

  push32(v, 0);    // 0: SEQ
  push8(v, 0);     // 4: Command
  push16(v, 1025); // 5: DDC port
  push16(v, 1026); // 7: DUC port
  push16(v, 1027); // 9: high priority from PC port
  push16(v, 1025); // 11: high priority to PC port
  push16(v, 1028); // 13: DDC audio port
  push16(v, 1029); // 15: DUC0 I&Q base port
  push16(v, 1035); // 17: DDC0 port (DDCx up to ... + 79)
  push16(v, 1026); // 19: Mic samples port
  push16(v, 1027); // 21: wideband ADC0 port
  push8(v, 0);     // 23: wideband enable
  push16(v, 512);  // 24: wideband samples per packet
  push8(v, 16); // 26: wideband sample size
  push8(v, 70); // 27: wideband update rate, ms per frame
  push8(v, 32); // 28: wideband packets per frame
  push16(v, 0);       // 29: memory mapped from PC port
  push16(v, 0);       // 31: memory mapped to PC port
  push16(v, 0);       // 33: reserved
  push16(v, 0);       // 35: reserved
  push8(v, 1 << 3); // 37: bit[3]=1 means frequency set as phase word
  push8(v, 1);  // 38: enable hardware reset / watchdog timer
  push8(v, 0);  // 39: big-endian, I&Q 3-byte format
  for(int i = 40; i < 56; i++)
    push8(v, 0);
  push8(v, 0);  // 56: Atlas/Mercury DDC config -- ignored?
  push8(v, 0);  // 57: ref clock source -- ignored?
  push8(v, tx_pa_);  // 58: enable PA
  int alexes = 1;
  if(do_both_)
    alexes |= 2; // enable alex1
  push8(v, alexes);  // 59: enable Alex filter boards, one bit each
  assert(v.size() == 60);

  assert(ntohs(radio_sin_.sin_port) == 1024);
  if(sendto(s_, v.data(), v.size(), 0, (sockaddr *) &radio_sin_,
            sizeof(radio_sin_)) < 0){
    perror("hpsdr: sendto");
  }
}

//
// send DDC setup packet to port 1025.
//
void
HPSDR::send_ddc()
{
  std::vector<char> v;

  push32(v, 0); // 0: seq
  push8(v, 2); // 4: two ADCs on the ANAN-100D
  push8(v, do_dither_?3:0); // 5: activate ADC dither, bitmap.
  push8(v, do_random_?3:0); // 6: activate ADC random, bitmap.

  unsigned int active = 0;
  for(int i = 0; i < NDDC && i < 8; i++){
    if(ddc_[i].active_)
      active |= (1 << i);
  }
  push8(v, active); // 7: activate DDCx, bitmap.

  for(int i = 8; i < 17; i++)
    push8(v, 0); // don't activate DDCs 8..79

  for(int unit = 0; unit < 80; unit++){
    int adc = (unit > 0 && do_both_) ? 1 : 0;
    push8(v, adc); // which ADC does DDCx talk to?
    int rate = unit < NDDC ? ddc_[unit].sdr_rate_ / 1000 : 48;
    push16(v, rate); // DDCx sampling rate; 48 is 48000
    push16(v, 0);  // reserved
    push8(v, 24); // DDCx I&Q sample size
  }
  for(int i = 497; i < 1443; i++)
    push8(v, 0);

  push8(v, 0); // 1443: not used
  assert(v.size() == 1444);

  sockaddr_in sin = radio_sin_;
  sin.sin_port = htons(1025);
  if(sendto(s_, v.data(), v.size(), 0, (sockaddr *) &sin,
            sizeof(sin)) < 0){
    perror("hpsdr: sendto");
  }
}

//
// send DUC setup packet to port 1026.
//
void
HPSDR::send_duc()
{
  std::vector<char> v;

  push32(v, 0); // 0: seq
  push8(v, 1); // 4: # of DACs
  push8(v, 0); // 5: 0 means not CW
  push8(v, 0); // 6: CW sidetone level
  push16(v, 0);// 7: CW sidetone hz
  push8(v, 0); // 9: CW keyer speed WPM
  push8(v, 50); // 10: CW weight
  push16(v, 0); // 11: CW hang time, ms
  push8(v, 0); // 13: RF delay, ms
  push16(v, 192); // 14: DUC sample rate
  push8(v, 24);  // 16: 16 bits/sample
  for(int i = 17; i <= 49; i++)
    push8(v, 0); // reserved
  push8(v, 0); // 50: mic control
  push8(v, 0); // 51: line in gain
  for(int i = 52; i <= 58; i++)
    push8(v, 0); // reserved
  push8(v, 0); // 59: reserved
  assert(v.size() == 60);

  sockaddr_in sin = radio_sin_;
  sin.sin_port = htons(1026);
  if(sendto(s_, v.data(), v.size(), 0, (sockaddr *) &sin,
            sizeof(sin)) < 0){
    perror("hpsdr: sendto");
  }
}

//
// compute 32-bit "phase word" corresponding to
// a specified frequency.
//
// phase_word[31:0] = 2^32 * frequency(Hz)/DSP clock frequency (Hz)
//
// DSP clock frequency is dependent on the Board type and is either
// specified in Appendix A or, if hardware specific, then as part of
// the Discovery response as specified in Appendix B.
//
unsigned int
HPSDR::ph(int hz)
{
  double x = hz / 122880000.0;
  x *= 4.0 * 1024.0 * 1024.0 * 1024.0; // 2^32
  return x;
}

//
// choose Alex pre-selector control bits
// based on min and max receive frequency,
// and on whether we're transmitting.
//
unsigned int
HPSDR::make_alex0()
{
  double lo = -1;
  double hi = -1;

  if(tx_active_){
    lo = tx_hz_;
    hi = tx_hz_;
  } else if(do_both_){
    // just DDC 0 on ADC 0
    hi = lo = ddc_[0].hz_;
  } else {
    // lowest and highest MHz we want to receive?
    // in case multiple DDCs.
    for(int unit = 0; unit < NDDC; unit++){
      int hz = ddc_[unit].hz_;
      if(hz > 0){
        if(lo < 0 || hz < lo){
          lo = hz;
        }
        if(hi < 0 || hz > hi){
          hi = hz;
        }
      }
    }
  }
  
  lo /= 1000000.0;
  hi /= 1000000.0;

  unsigned int alex0 = 0;

  //
  // TX side.
  //

  alex0 |= 1 << 24; // Ant 1

  if(do_lpf_){
    // low-pass filter (shared with TX).
    if(hi < 2.0){
      alex0 |= 1 << 23; // 160M
    } else if(hi < 4.0){
      alex0 |= 1 << 22; // 80M
    } else if(hi < 7.3){
      alex0 |= 1 << 21; // 60/40M
    } else if(hi < 14.350){
      alex0 |= 1 << 20; // 30/20M
    } else if(hi < 21.450){
      alex0 |= 1 << 31; // 17/15M
    } else if(hi < 30.0){
      alex0 |= 1 << 30; // 12/10M
    } else {
      alex0 |= 1 << 29; // 6M/bypass
    }
  } else {
    alex0 |= 1 << 29; // 6M/Bypass LPF
  }

  //
  // RX side.
  //

  if(tx_active_){
    // do nothing -- T/R switch will disconnect the receivers.
  } else {
    if(do_ext1_){
      alex0 |= 1 << 9; // Ext 1
      if(board_type_ == ANGELIA_BOARD){
        alex0 |= 1 << 11; // ByPass (disconnect Ant 1 from RX)
      }
      if(board_type_ == ORION2_BOARD){
        alex0 |= 1 << 14; // Rx Master in Sel.
      }
    }

    if(board_type_ == ANGELIA_BOARD){
      int att = do_attenuate_;
      if(att >= 20){
        alex0 |= 1 << 13; // 20 dB Alex attenuator.
        att -= 20;
      }
      if(att >= 10){
        alex0 |= 1 << 14; // 10 dB Alex attenuator.
        att -= 10;
      }
    }

    if(board_type_ == ANGELIA_BOARD){
      if(do_hpf_){
        // high-pass filter.
        if(lo >= 50.0){
          alex0 |= 1 << 3; // 6M with LNA
        } else if(lo >= 20.0){
          alex0 |= 1 << 2;
        } else if(lo >= 13){
          alex0 |= 1 << 1;
        } else if(lo >= 9.5){
          alex0 |= 1 << 4;
        } else if(lo >= 6.5){
          alex0 |= 1 << 5;
        } else if(lo >= 1.5){
          alex0 |= 1 << 6;
        } else {
          alex0 |= 1 << 12; // HF Bypass (maybe doesn't work?)
        }
      } else {
        alex0 |= 1 << 6; // 1.5 MHz HPF
      }
    }
    if(board_type_ == ORION2_BOARD){
      if(do_hpf_){
        // band-pass filters
        if(lo >= 1.5 && hi < 2.1){
          alex0 |= 1 << 6;
        } else if(lo >= 2.1 && hi < 5.5){
          alex0 |= 1 << 5;
        } else if(lo >= 5.5 && hi < 10.99){
          alex0 |= 1 << 4;
        } else if(lo >= 11.0 && hi < 21){
          alex0 |= 1 << 1;
        } else if(lo >= 21 && hi < 27){
          alex0 |= 1 << 2;
        } else if(lo >= 27 && hi < 61.4){
          alex0 |= 1 << 3;
        } else {
          alex0 |= 1 << 12; // Bypass
        }
      } else {
        alex0 |= 1 << 12; // Bypass
      }
    }
  }

  return alex0;
}

//
// 2nd Alex board, connected to ADC 1.
// only used if do_both_.
// and in that case, units 1, 2, and 3
// use ADC1 and Alex1.
// Alex1 has only RX bandpass filters (no TX filters).
//
unsigned int
HPSDR::make_alex1()
{

  if(board_type_ != ORION2_BOARD){
    // ANAN-100D has no Alex 1.
    return 0;
  }

  if(tx_active_ || do_both_ == 0){
    // RX2 ground
    return 1 << 8;
  }
  
  // lowest and highest MHz we want to receive?
  // in case multiple DDCs.
  double lo = -1;
  double hi = -1;
  for(int unit = 1; unit < NDDC; unit++){
    int hz = ddc_[unit].hz_;
    if(hz > 0){
      if(lo < 0 || hz < lo){
        lo = hz;
      }
      if(hi < 0 || hz > hi){
        hi = hz;
      }
    }
  }
  
  lo /= 1000000.0;
  hi /= 1000000.0;

  unsigned int alex1 = 0;

  if(do_hpf_){
    // band-pass filters
    if(lo >= 1.5 && hi < 2.1){
      alex1 |= 1 << 6;
    } else if(lo >= 2.1 && hi < 5.5){
      alex1 |= 1 << 5;
    } else if(lo >= 5.5 && hi < 10.99){
      alex1 |= 1 << 4;
    } else if(lo >= 11.0 && hi < 21){
      alex1 |= 1 << 1;
    } else if(lo >= 21 && hi < 27){
      alex1 |= 1 << 2;
    } else if(lo >= 27 && hi < 61.4){
      alex1 |= 1 << 3;
    } else {
      alex1 |= 1 << 12; // Bypass
    }
  } else {
    alex1 |= 1 << 12; // Bypass
  }

  return alex1;
}

//
// band selector bits on OC0..OC3 on the DB9 connector
// on the back of the ANAN-7000dle.
//
unsigned int
HPSDR::make_band()
{
  int bits = 0;
  double mhz = tx_hz_ / 1000000.0;

  if(mhz < 1.8){
    // ???
  } else if(mhz < 2.0){
    bits = 1; // 160M
  } else if(mhz < 4.0){
    bits = 2; // 80M
  } else if(mhz < 7.3){
    bits = 3; // 40M
  } else if(mhz < 10.2){
    bits = 4; // 30M
  } else if(mhz < 14.350){
    bits = 5; // 20M
  } else if(mhz < 19){
    bits = 6; // 17M
  } else if(mhz < 22){
    bits = 7; // 15M
  } else if(mhz < 25){
    bits = 8; // 12M
  } else if(mhz < 30){
    bits = 9; // 10M
  } else if(mhz < 55){
    bits = 10; // 6M
  }

  bits = 3; // XXX

  bits <<= 1;

  return bits;
}

//
// give a desired PA power output, return a good drive
// level to put in the high-priority packets.
//
// I don't know how to do this correctly. All I know is
// the drive level required to get about 25 watts on
// my ANAN-100D.
//
// the drive-vs-power relationship is not linear,
// but I don't know what it is.
//
int
HPSDR::watts2drive(int hz, int watts)
{
  if(watts > 30)
    watts = 30;
  int drive = 10; // pretty low.
  double mhz = hz / 1000000.0;

  if(board_type_ == ORION2_BOARD){
    // use Blue ANAN-7000DLE measurements at 5 watts
    if(mhz >= 1.5 && mhz <= 2){
      drive = 23 * watts / 5.0;
    } else if(mhz >= 3.5 && mhz <= 4){
      drive = 21 * watts / 5.0;
    } else if(mhz >= 7 && mhz <= 8){
      drive = 25 * watts / 5.0;
    } else if(mhz >= 10 && mhz <= 11){
      drive = 26 * watts / 5.0;
    } else if(mhz >= 14 && mhz <= 15){
      drive = 30 * watts / 5.0;
    } else if(mhz >= 18 && mhz <= 19){
      drive = 24 * watts / 5.0;
    } else if(mhz >= 21 && mhz <= 22){
      drive = 27 * watts / 5.0;
    } else if(mhz >= 24 && mhz <= 25){
      drive = 30 * watts / 5.0;
    } else if(mhz >= 28 && mhz <= 30){
      drive = 37 * watts / 5.0;
    } else if(mhz >= 50 && mhz <= 54){
      drive = 37 * watts / 5.0;
    } else {
      static int once = 0;
      if(once == 0){
        fprintf(stderr, "HPSDR::watts2drive(%d, %d) : %d hz not recognized\n",
                hz, watts, hz);
        once = 1;
      }
    }
  }

  if(board_type_ == ANGELIA_BOARD){
    if(watts <= 10){
      // use ANAN-100D measurements at 5 watts
      if(mhz >= 1.5 && mhz <= 2){
        drive = 16 * watts / 5.0;
      } else if(mhz >= 3.5 && mhz <= 4){
        drive = 17 * watts / 5.0;
      } else if(mhz >= 7 && mhz <= 8){
        drive = 20 * watts / 5.0;
      } else if(mhz >= 10 && mhz <= 11){
        drive = 21 * watts / 5.0;
      } else if(mhz >= 14 && mhz <= 15){
        drive = 26 * watts / 5.0;
      } else if(mhz >= 18 && mhz <= 19){
        drive = 26 * watts / 5.0;
      } else if(mhz >= 21 && mhz <= 22){
        drive = 30 * watts / 5.0;
      } else if(mhz >= 24 && mhz <= 25){
        drive = 25 * watts / 5.0;
      } else if(mhz >= 28 && mhz <= 30){
        drive = 31 * watts / 5.0;
      } else if(mhz >= 50 && mhz <= 54){
        drive = 55 * watts / 5.0;
      } else {
        static int once = 0;
        if(once == 0){
          fprintf(stderr, "HPSDR::watts2drive(%d, %d) : %d hz not recognized\n",
                  hz, watts, hz);
          once = 1;
        }
      }
    } else {
      // use ANAN-100D measurements at 25 watts
      if(mhz >= 1.5 && mhz <= 2){
        drive = 27 * watts / 25.0;
      } else if(mhz >= 3.5 && mhz <= 4){
        drive = 30 * watts / 25.0;
      } else if(mhz >= 7 && mhz <= 8){
        drive = 33 * watts / 25.0;
      } else if(mhz >= 10 && mhz <= 11){
        drive = 35 * watts / 25.0;
      } else if(mhz >= 14 && mhz <= 15){
        drive = 45 * watts / 25.0;
      } else if(mhz >= 18 && mhz <= 19){
        drive = 42 * watts / 25.0;
      } else if(mhz >= 21 && mhz <= 22){
        drive = 50 * watts / 25;
      } else if(mhz >= 24 && mhz <= 25){
        drive = 40 * watts / 25;
      } else if(mhz >= 28 && mhz <= 30){
        drive = 50 * watts / 25;
      } else if(mhz >= 50 && mhz <= 54){
        drive = 100 * watts / 25;
      } else {
        static int once = 0;
        if(once == 0){
          fprintf(stderr, "HPSDR::watts2drive(%d, %d) : %d hz not recognized\n",
                  hz, watts, hz);
          once = 1;
        }
      }
    }
  }

  return drive;
}

//
// send high-priority setup packet to port 1027.
//
void
HPSDR::send_high()
{
  std::vector<char> v;

  push32(v, 0); // 0: seq
  push8(v, run_ | (tx_ptt_ << 1));  // 4: bit 0 = enable SDR (run), bit 1 = enable TX
  push8(v, 0);  // 5: CW
  push8(v, 0);  // 6: reserved
  push8(v, 0);  // 7: reserved
  push8(v, 0);  // 8: reserved
  for(int i = 0; i < 80; i++){
    // 9: DDC0 frequency/phase word
    if(i < NDDC && ddc_[i].hz_){
      push32(v, ph(ddc_[i].hz_));
    } else {
      push32(v, 0);
    }
  }
  push32(v, ph(tx_hz_)); // 329: DUC0 frequency/phase word
  push32(v, 0); // 333: DUC1 frequency/phase word
  push32(v, 0); // 337: DUC2 frequency/phase word
  push32(v, 0); // 341: DUC3 frequency/phase word
  if(tx_drive_ > 0 && tx_drive_override_ >= 0){ 
    push8(v, tx_drive_override_); // 345: DUC0 drive level
  } else {
    push8(v, tx_drive_); // 345: DUC0 drive level
  }
  push8(v, 0); // 346: DUC1 drive level
  push8(v, 0); // 347: DUC2 drive level
  push8(v, 0); // 348: DUC3 drive level
  for(int i = 349; i < 1400; i++)
    push8(v, 0); // ???
  push8(v, 0); // 1400: transverter and audio enable
  push8(v, make_band()); // 1401: open collector outputs
  push8(v, 0); // 1402: user outputs DB9 pins 1-4
  push8(v, 0); // 1403: Mercury attenuator
  for(int i = 1404; i < 1430; i++)
    push8(v, 0); // reserved
  unsigned int alex1 = make_alex1();
  push16(v, alex1); // 1430: Alex 1
  unsigned int alex0 = make_alex0();
  push32(v, alex0); // 1432: Alex 0
  for(int i = 1436; i <= 1441; i++)
    push8(v, 0);
  push8(v, step_[1]); // 1442: step attenuator 1
  push8(v, step_[0]); // 1443: step attenuator 0
  assert(v.size() == 1444);

  sockaddr_in sin = radio_sin_;
  sin.sin_port = htons(1027);
  if(sendto(s_, v.data(), v.size(), 0, (sockaddr *) &sin,
            sizeof(sin)) < 0){
    perror("hpsdr: sendto");
  }
}

//
// run one I/Q sample through the low-pass filter.
// (preparatory to down-sampling from sdr_rate_ to rate_.)
void
HPSDR::lowpass(int unit, double &ii, double &qq)
{
  liquid_float_complex x, y;
  x.real = ii;
  x.imag = qq;
  firfilt_crcf_push(ddc_[unit].filter_, x);
  firfilt_crcf_execute(ddc_[unit].filter_, &y);
  ii = y.real;
  qq = y.imag;
}

void
HPSDR::process_data(const std::vector<char> &buf, int unit)
{
  assert(unit >= 0 && unit < NDDC);

  if(ddc_[unit].active_ == 0)
    return;

  // 4..12 is 64-bit timestamp

  int bits_per_sample = r16(buf, 12); // 24
  if(bits_per_sample != 24){
    static int warned = 0;
    if(warned == 0){
      warned = 1;
      fprintf(stderr, "HPSDR::process_data(): unexpected %d bits/sample\n",
              bits_per_sample);
      return;
    }
  }

  int samples_per_frame = r16(buf, 14); // typically 238
  if(samples_per_frame * 6 + 16 > (int) buf.size()){
    fprintf(stderr, "oops samples_per_frame %d buf.size() %d\n",
            samples_per_frame, (int) buf.size());
    return;
  }

  std::vector<std::complex<double>> o;

  unsigned int seq = r32(buf, 0);
  if(seq != ddc_[unit].expect_seq_){
    fprintf(stderr, "HPSDR unit %d got seq %d expected %d\n",
            unit, seq, ddc_[unit].expect_seq_);
    // fake the missing samples.
    int needed = (ddc_[unit].expect_seq_ - seq) * samples_per_frame;
    for(int i = 0; i < needed; i++){
      if((ddc_[unit].count_ % (ddc_[unit].sdr_rate_ / ddc_[unit].rate_)) == 0){
        o.push_back(std::complex<double>(0, 0));
      }
      ddc_[unit].count_ += 1;
    }
  }
  ddc_[unit].expect_seq_ = seq + 1;

  // 16: now 3 bytes Q, 3 bytes I, ...

  for(int i = 0; i < samples_per_frame; i++){
    int off = i * 6 + 16;
    // Q comes first, then I
    double qq = r24(buf, off + 0); // signed...
    double ii = r24(buf, off + 3); // signed...

    if(ddc_[unit].rate_ < ddc_[unit].sdr_rate_){
      // low-pass filter, preparatory to rate reduction.
      lowpass(unit, ii, qq);
    }

    if((ddc_[unit].count_ % (ddc_[unit].sdr_rate_ / ddc_[unit].rate_)) == 0){
      o.push_back(std::complex<double>(ii, qq));
    }

    ddc_[unit].count_ += 1;
  }

  // mu_ protects input_ from concurrent get()
  ddc_[unit].mu_.lock();
  ddc_[unit].input_.insert(ddc_[unit].input_.end(), o.begin(), o.end());
  ddc_[unit].mu_.unlock();
}

void
HPSDR::process_priority_status(const std::vector<char> &buf,
                               const sockaddr_in &from)
{
  double tt = now();

  for(int adc = 0; adc < 2; adc++){
    if(tx_active_ == 0){
      if(buf[5] & (1 << adc)){
        if(tt > last_clipped_[adc] + 5){
          last_clipped_[adc] = tt;
          last_unclipped_[adc] = tt;
          step_[adc] = std::min(step_[adc] + 10, 31);
          fprintf(stderr, "hpsdr: clip %d -> %d dB\n", adc, step_[adc]);
        }
      } else if(tt > last_unclipped_[adc] + 60){
        step_[adc] = std::max(0, step_[adc] - 5);
        last_unclipped_[adc] = tt;
      }
    }
  }

  // voltage and power conversions from Appendix A.

  double vconstant = (board_type_ == ORION2_BOARD ? 5.0 : 3.3);
  double constant1 = (board_type_ == ORION2_BOARD ? 5.0 : 3.3);
  double constant2 = (board_type_ == ORION2_BOARD ? 0.108 : 0.095);

  int raw_exciter_power = r16(buf, 6);
  double exciter_power = pow((raw_exciter_power / 4095.0) * constant1, 2)
    / constant2;

  int raw_forward_power = r16(buf, 14);
  double forward_power = pow((raw_forward_power / 4095.0) * constant1, 2)
    / constant2;

  int raw_reverse_power = r16(buf, 22);
  double reverse_power = pow((raw_reverse_power / 4095.0) * constant1, 2)
    / constant2;

  double swr = 0.0;
  if(forward_power > 0){
    double a = 1 + sqrt(reverse_power / forward_power);
    double b = 1 - sqrt(reverse_power / forward_power);
    if(b != 0){
      swr = a / b;
    }
  }

  int raw_supply_voltage = r16(buf, 49);
  double supply_voltage = (raw_supply_voltage / 4095.0) * vconstant;

  if(tx_active_ && tt > last_swr_print_ + 4){
    fprintf(stderr, "exciter %.1f mW, fwd %.1f W, rev %.1f W, swr %.1f, supply %.1f\n",
            exciter_power,
            forward_power,
            reverse_power,
            swr,
            supply_voltage);
    last_swr_print_ = tt;
  }

  // the 3.8 document says these two are no longer supported,
  // which is too bad.
  if(buf[4] & (1 << 5)){
    fprintf(stderr, "FIFO almost empty\n");
  }
  if(buf[4] & (1 << 6)){
    fprintf(stderr, "FIFO almost full\n");
  }
}
  
void
HPSDR::process_packet(const std::vector<char> &buf,
                      const sockaddr_in &from)
{
  int sport = ntohs(from.sin_port);

  if(state_ == Discovering && buf.size() == 60 &&
     (buf[4] == 0x02 || buf[4] == 0x03)){
    // discovery response
    // 0x02 means radio is idle
    // 0x03 means radio is already talking to someone
    process_discovery_response(buf, from);
    send_general();
    send_ddc();
    send_duc();
    send_high();
    last_general_ = now();
    state_ = Running;
  } else if(state_ == Running && sport == 1025){
    // high priority status
    process_priority_status(buf, from);
  } else if(state_ == Running && sport == 1026){
    // mic samples
  } else if(state_ == Running && sport == 1027){
    // wideband samples
  } else if(state_ == Running && (sport >= 1035 && sport < 1035+NDDC)){
    // DDC0 I&Q samples
    process_data(buf, sport - 1035);
  } else {
    fprintf(stderr, "unexpected packet from %s:%d\n",
            inet_ntoa(from.sin_addr),
            sport);
  }
}

//
// fetch all the buffered input.
// does not block.
//
void
HPSDR::get(int unit, std::vector<std::complex<double>> &buf)
{
  assert(unit >= 0 && unit < NDDC && ddc_[unit].active_);
  buf.clear();
  ddc_[unit].mu_.lock(); // lock against concurrent append by loop().
  buf.swap(ddc_[unit].input_);
  ddc_[unit].mu_.unlock();
}

//
// send one packet of 240 I/Q samples.
// 192000 samples/second, 24 bits.
// expects iq.real() and iq.imag() to be -1 .. 1,
// scales to 24-bit fixnum.
//
void
HPSDR::send_to_duc(const std::vector<std::complex<double>> &iq, int off)
{
  std::vector<char> v;

  push32(v, tx_seq_++);
  for(int i = 0; i < 240; i++){
    if(i + off < (int) iq.size()){
      double x = iq[i+off].imag();
      if(x < -1) x = -1;
      if(x > 1) x = 1;
      x *= 0.99;
      x *= 1 << 23;
      push24(v, x);

      x = iq[i+off].real();
      if(x < -1) x = -1;
      if(x > 1) x = 1;
      x *= 0.99;
      x *= 1 << 23;
      push24(v, x);
    }
  }

  assert(v.size() == 1444);

  sockaddr_in sin = radio_sin_;
  sin.sin_port = htons(1029); // DAC0
  if(sendto(s_, v.data(), v.size(), 0, (sockaddr *) &sin,
            sizeof(sin)) < 0){
    perror("hpsdr: sendto");
  }
}

//
// transmit.
// hz is DUC frequency.
// v[] should be I/Q at rate samples/second, each -1..1
// sets PTT, sends samples to DUC, clears PTT.
//
// set DUC frequency -- done
// set proper low-pass filter -- done
// DUC drive level in high priority packet -- done
// PTT bit in high priority packet -- done
// enable PA in byte 58 of general packet -- not done
// do I need to explicitly set T/R control (Alex bit 27)? no
// watch full/empty status in High Priority Status packet?
//
void
HPSDR::tx_iq(int hz, const std::vector<std::complex<double>> &iq, int rate)
{
  assert(rate <= 192000 && rate > 0);

  fprintf(stderr, "tx_iq hz=%d n=%d rate=%d power=%d drive=%d pa=%d\n",
          hz, (int) iq.size(), rate, tx_power_,
          watts2drive(hz, tx_power_), tx_pa_);

  // step 1:
  // set DUC frequency.
  // set up Alex relays for TX and mute RX.
  // (no PTT, no drive)
  // delay.
  tx_active_ = 1;
  tx_ptt_ = 0;
  tx_hz_ = hz;
  tx_drive_ = 0;
  send_high();
  usleep(10 * 1000);

  // step 2:
  // set drive level.
  // set PTT.
  // delay again.
  tx_ptt_ = 1;
  send_high();
  usleep(10 * 1000);

  // step 3:
  // non-zero drive.
  tx_drive_ = watts2drive(tx_hz_, tx_power_);
  send_high();

  // what does PTT do? PA? T/R switch?
  // do I need to delay after setting PTT?
  // when should I set the drive level?
  // when should I start sending samples?
  // what happens if ptt and drive are set but there are
  //   no I/Q samples in the DUC FIFO?
  // should something ramp up?

  double t0 = now();
  std::vector<std::complex<double>> buf;
  int oi = 0; // output index
  while(1){
    double iir = (oi / 192000.0) * rate;
    int ii = iir; // truncate; iir >= ii
    if(ii >= (int) iq.size())
      break;
    std::complex<double> c;
    if(ii >= (int) iq.size() - 1){
      c = iq[iq.size() - 1];
    } else {
      double f = iir - ii;
      c = (iq[ii] * (1.0 - f)) + (iq[ii+1] * f);
    }

    buf.push_back(c);
    oi += 1;

    if(buf.size() >= 240){
      send_to_duc(buf, 0);
      buf.clear();

      // the FIFO feeding the DUC probably has 4096 I/Q-bit entries.
      // sleep until there are likely to be only 2048 entries left in the FIFO.
      if(oi >= 2048){
        double t1 = now();
        double tt = t0 + (oi - 2048) / 192000.0;
        if(t1 < tt){
          usleep((tt - t1) * 1000000);
        }
      }
    }
  }

  // try to wait until the radio has consumed and transmitted everything.
  double t1 = now();
  double tt = t0 + oi / 192000.0;
  if(t1 < tt){
    usleep((tt - t1) * 1000000);
  }

  // termination stage 1:
  // zero drive.
  tx_drive_ = 0;
  send_high();
  usleep(10 * 1000);

  // stage 2:
  // turn off PTT.
  tx_ptt_ = 0;
  send_high();
  usleep(10 * 1000);

  // stage 3:
  // enable the receiver.
  tx_active_ = 0;
  send_high();
}

//
// call with ordinary audio samples at rate.
// samples should be -1..1.
// DUC center at hz.
//
void
HPSDR::tx_real(int hz, const std::vector<double> &v, int rate)
{
  // FFTW_MEASURE is too expensive to call when transmitting!
  extern int fftw_type;
  int otype = fftw_type;
  fftw_type = FFTW_ESTIMATE;

  std::vector<std::complex<double>> vv = analytic(v, "tx_real");
  std::vector<double> q = vimag(vv);

  std::vector<std::complex<double>> iq(q.size());
  for(int i = 0; i < (int) v.size(); i++){
    iq[i] = std::complex(v[i], q[i]);
  }
  tx_iq(hz, iq, rate);

  fftw_type = otype;
}

//
// call this in its own thread.
//
void
HPSDR::loop()
{
  while(1){
    tick();
    assert(state_ != None);
    if(wait_socket(s_, 200))
      process();
  }
}
