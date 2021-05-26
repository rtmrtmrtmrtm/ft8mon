//
// RFSpace SDR-IP, NetSDR, CloudIQ, and CloudSDR, in I/Q mode.
//
// Robert Morris, AB1HL.
//

#include "sdrip.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <string>
#include "util.h"

SDRIP::SDRIP(int rate)
{
  rate_ = rate;
  tcp_ = -1;
  udp_ = -1;
  adc_rate_ = -1;
  control_th_ = 0;
  data_th_ = 0;
  iqcount_ = 0;
  outcount_ = 0;
}


void
SDRIP::open(std::string chan)
{
  std::string hostname;
  
  long long hz = -1;
  int comma = chan.find(",");
  if(comma >= 0){
    hz = atof(chan.c_str() + comma + 1) * 1000000.0;
    hostname = std::string(chan, 0, comma);
  } else {
    hostname = chan;
  }

  int port = 50000;
  int colon = hostname.find(":");
  if(colon >= 0){
    port = atoi(hostname.c_str() + colon + 1);
    hostname = std::string(hostname, 0, colon);
  }

  tcp_ = socket(AF_INET, SOCK_STREAM, 0);
  assert(tcp_ >= 0);

  sockaddr_in sin;
  memset(&sin, 0, sizeof(sin));
  sin.sin_family = AF_INET;
  sin.sin_port = htons(port);
  sin.sin_addr.s_addr = inet_addr(hostname.c_str());

  int ret = connect(tcp_, (sockaddr *) &sin, sizeof(sin));
  if(ret < 0){
    fprintf(stderr, "SDRIP: connect(%s) failed\n", hostname.c_str());
    exit(1);
  }

  // ask the kernel to allocate us a UDP port number on which
  // to listen to incoming I/Q data.
  udp_ = socket(AF_INET, SOCK_DGRAM, 0);
  assert(udp_ >= 0);
  sockaddr_in usin;
  memset(&usin, 0, sizeof(usin));
  usin.sin_family = AF_INET;
  if(bind(udp_, (sockaddr *) &usin, sizeof(usin)) < 0){
    perror("SDRIP: udp bind");
    exit(1);
  }
  socklen_t usin_len = sizeof(usin);
  if(getsockname(udp_, (sockaddr *) &usin, &usin_len) < 0){
    perror("SDRIP: getsockname");
    exit(1);
  }
  sockaddr_in tsin;
  socklen_t tsin_len = sizeof(tsin);
  if(getsockname(tcp_, (sockaddr *) &tsin, &tsin_len) < 0){
    perror("SDRIP: getsockname");
    exit(1);
  }
  int udp_port = ntohs(usin.sin_port);
  unsigned int udp_ip = ntohl(tsin.sin_addr.s_addr);

  control_th_ = new std::thread( [ this ] () { control_reader(); } );
  data_th_ = new std::thread( [ this ] () { data_reader(); } );

  // if we send requests too soon, the SDRIP ignores them.
  usleep(100*1000);

  set_state(1); // stop
  adc_rate_ = get_adc_rate();
  if(hz >= 0){
    set_frequency(hz);
  }
  if(is_sdrip()){
    iq_rate_ = 32000;
  } else if(is_cloudsdr()){
    if((122880000 % (rate_ * 4)) == 0){
      // CloudSDR can directly generate rate_.
      iq_rate_ = rate_;
    } else if(rate_ <= 32000){
      iq_rate_ = 32000;
    } else {
      // should use one of the radio's higher sample rates.
      assert(0);
    }
  } else {
    assert(0);
  }
  set_iq_rate(iq_rate_);
  set_udp_addr(udp_ip, udp_port);
  set_gain(0);
  set_packet_size(0);
  // set_ad_mode(1);
  set_state(2); // run

  if(rate_ < iq_rate_){
    // Liquid DSP anti-alias filter for rate conversion.
    int h_len = estimate_req_filter_len(0.01, 60.0);
    float h[h_len];
    double cutoff = (rate_ / (double) iq_rate_) / 2.0;
    
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
    
    filter_ = firfilt_crcf_create(h, h_len);
  }

  printf("name: %s\n", get_name().c_str());
  printf("serial: %s\n", get_serial().c_str());
  printf("status: %s\n", get_status().c_str());
  printf("frequency: %lld\n", get_frequency());
  printf("gain: %d\n", get_gain());
  printf("A/D mode: %d\n", get_ad_mode());
  printf("rf filter: %d\n", get_rf_filter());
  printf("packet size: %d\n", get_packet_size());
  printf("ADC rate: %d\n", get_adc_rate());
  printf("I/Q rate: %d\n", get_iq_rate());
  {
    unsigned int ip;
    int port;
    get_udp_addr(ip, port);
    printf("UDP Addr: %08x %d\n", ip, port);
  }
  fflush(stdout);
}

std::string
SDRIP::get_name()
{
  std::vector<int> v = getitem(0x0001);
  std::string s;
  for(int i = 4; i < (int) v.size(); i++){
    if(v[i] != '\0')
      s.push_back(v[i]);
  }
  return s;
}

std::string
SDRIP::get_serial()
{
  std::vector<int> v = getitem(0x0002);
  std::string s;
  for(int i = 4; i < (int) v.size(); i++){
    if(v[i] != '\0')
      s.push_back(v[i]);
  }
  return s;
}

std::string
SDRIP::get_status()
{
  std::vector<int> v = getitem(0x0005);
  std::string s;
  for(int i = 4; i < (int) v.size(); i++){
    if(v[i] == 0x0B){
      s = s + "Idle ";
    } else if(v[i] == 0x0C){
      s = s + "Busy ";
    } else if(v[i] == 0x0E){
      s = s + "BootIdle ";
    } else if(v[i] == 0x0F){
      s = s + "BootBusy ";
    } else if(v[i] == 0x20){
      s = s + "Overload ";
    } else if(v[i] == 0x80){
      s = s + "Error ";
    } else {
      s = s + "?? ";
    }
  }
  return s;
}

//
// hz
//
long long
SDRIP::get_frequency()
{
  std::vector<int> v = getitem(0x0020, 0);

  if(v.size() < 10)
    return -1;
  
  long long hz = 0;
  // LSB first
  for(int i = 0; i < 5; i++){
    hz |= ((v[i+5] & 0xff) << (8 * i));
  }
  
  return hz;
}

void
SDRIP::set_frequency(long long hz)
{
  std::vector<int> data(6);
  data[0] = 0; // channel ID
  for(int i = 0; i < 5; i++){
    // least significant byte first
    data[i+1] = hz & 0xff;
    hz >>= 8;
  }
  setitem(0x0020, data);
}

// 1 = stop
// 2 = run
void
SDRIP::set_state(int state)
{
  std::vector<int> data;
  data.push_back(0x80); // I/Q data
  data.push_back(state);
  data.push_back(0x80); // 24-bit data, contiguous
  data.push_back(0x00);
  setitem(0x0018, data);
}

void
SDRIP::set_iq_rate(int rate)
{
  if(is_cloudsdr()){
    assert((122880000 % rate) == 0);
    assert((122880000 % (rate * 4)) == 0);
  }
  if(is_sdrip()){
    assert(rate == 32000);
  }
  std::vector<int> data;
  data.push_back(0x00); // ignored
  data.push_back((rate >> 0) & 0xff);
  data.push_back((rate >> 8) & 0xff);
  data.push_back((rate >> 16) & 0xff);
  data.push_back((rate >> 24) & 0xff);
  setitem(0x00B8, data);
}

int
SDRIP::get_iq_rate()
{
  std::vector<int> v = getitem(0x00B8, 0);
  assert(v.size() == 9);
  int r = 0;
  r |= (v[5] << 0);
  r |= (v[6] << 8);
  r |= (v[7] << 16);
  r |= (v[8] << 24);
  return r;
}

//
// tell the SDRIP where to send UDP I/Q data packets.
// ip and port are host byte order.
//
void
SDRIP::set_udp_addr(unsigned int ip, int port)
{
  std::vector<int> data;
  data.push_back((ip >> 0) & 0xff);
  data.push_back((ip >> 8) & 0xff);
  data.push_back((ip >> 16) & 0xff);
  data.push_back((ip >> 24) & 0xff);
  data.push_back((port >> 0) & 0xff);
  data.push_back((port >> 8) & 0xff);
  setitem(0x00C5, data);
}

void
SDRIP::get_udp_addr(unsigned int &ip, int &port)
{
  std::vector<int> v = getitem(0x00C5, -1);
  ip = 0;
  ip |= (v[4] << 0);
  ip |= (v[5] << 8);
  ip |= (v[6] << 16);
  ip |= (v[7] << 24);
  port = 0;
  port |= (v[8] << 0);
  port |= (v[9] << 8);
}

//
// 0, -10, -20, -30
//
int
SDRIP::get_gain()
{
  std::vector<int> v = getitem(0x0038, 0);
  if(v.size() < 6)
    return -1;
  return (signed char) v[5];
}

void
SDRIP::set_gain(int gain)
{
  std::vector<int> data;
  data.push_back(0);
  data.push_back(gain);
  setitem(0x0038, data);
}

// 2 = A/D gain 1.5, 0 = A/D gain 1.0
int
SDRIP::get_ad_mode()
{
  std::vector<int> v = getitem(0x008A, 0);
  assert(v.size() == 6);
  return v[5];
}

void
SDRIP::set_ad_mode(int mode)
{
  std::vector<int> data;
  data.push_back(0);
  data.push_back(mode);
  setitem(0x008A, data);
}

int
SDRIP::get_rf_filter()
{
  std::vector<int> v = getitem(0x0044, 0);
  if(v.size() < 6)
    return -1;
  return v[5] & 0xff;
}

//
// roughly 122880000 for CloudIQ/CloudSDR, 80000000 for SDRIP/NetSDR
//
int
SDRIP::get_adc_rate()
{
  std::vector<int> v = getitem(0x00B0, 0);
  assert(v.size() == 9);
  int rate = 0;
  rate |= (v[5] << 0);
  rate |= (v[6] << 8);
  rate |= (v[7] << 16);
  rate |= (v[8] << 24);
  return rate;
}

// 0 = large UDP packets
// 1 = small UDP packets
int
SDRIP::get_packet_size()
{
  std::vector<int> v = getitem(0x00C4, -1);
  assert(v.size() == 5);
  return v[4];
}

void
SDRIP::set_packet_size(int sz)
{
  std::vector<int> data;
  data.push_back(sz);
  setitem(0x00C4, data);
}

//
// watch the TCP connection for replies to control messages.
//
void
SDRIP::control_reader()
{
  while(1){
    std::vector<int> v = readone();

    int type = (v[1] >> 5) & 0x7;
    if(type == 0){
      // response to set or request item
      reply_mu_.lock();
      reply_ = v;
      reply_mu_.unlock();
    } else if(v.size() == 7 && type == 4){
      // a seven-byte message starting like this is squelched data:
      // 07 80 ...
      // such data starts arriving as soon as a connection is made.
    } else {
      fprintf(stderr, "SDRIP unknown type=%d, ", type);
      for(int i = 0; i < 6 && i < (int) v.size(); i++){
        fprintf(stderr, "%02x ", v[i]);
      }
      fprintf(stderr, "\n");
    }
  }
}

//
// process one I/Q sample from a UDP packet.
//
void
SDRIP::one_iq(const std::complex<double> iq0)
{
  if(rate_ < iq_rate_){
    // anti-alias filter
    liquid_float_complex x, y;
    x.real = iq0.real();
    x.imag = iq0.imag();
    firfilt_crcf_push(filter_, x);
    firfilt_crcf_execute(filter_, &y);
    std::complex<double> iq1(y.real, y.imag);
  
    // convert iq_rate_ samples/second to rate_.
    iqcount_ += 1;
    if((long long)((iqcount_ / (double) iq_rate_) * rate_) > outcount_){
      samples_.push_back(iq1);
      outcount_ += 1;
    }
  } else {
    samples_.push_back(iq0);
  }
}

//
// read I/Q UDP data packets from the SDRIP, append
// the data to samples_.
//
void
SDRIP::data_reader()
{
  int next_seq = 0;
  unsigned char buf[2048];
  while(1){
    int cc = recv(udp_, (char*)buf, sizeof(buf), 0);
    if(cc < 0){
      perror("sdrip: recv");
      exit(1);
    }
    // 0xA4 0x85 -- 1440 data bytes, 240 24-bit I/Q samples
    // 0x84 0x81 -- 384 data bytes, 64 24-bit I/Q samples
    if(cc < 4 || (buf[0] != 0xA4 && buf[0] != 0x84) || (buf[1] != 0x85 && buf[1] != 0x81)){
      fprintf(stderr, "sdrip weird UDP packet\n");
      continue;
    }

    int nmissing = 0;
    int seq = buf[2] | (buf[3] << 8);
    if(seq != next_seq){
      fprintf(stderr, "sdrip: expected seq %d, got %d\n", next_seq, seq);
      if(seq > next_seq){
        nmissing = 240 * (seq - next_seq);
      }
    }
    if(seq == 0xffff){
      // wrap around to 1, not 0
      next_seq = 1;
    } else {
      next_seq = seq + 1;
    }

    samples_mu_.lock();

    for(int i = 0; i < nmissing; i++){
      std::complex<double> junk(0, 0);
      one_iq(junk);
    }

    for(int i = 4; i < cc; i += 6){
      int ii = buf[i+0] | (buf[i+1] << 8) | (buf[i+2] << 16);
      if(ii & (1 << 23)){
        // sign-extend
        ii |= 0xff000000;
      }
      int qq = buf[i+3] | (buf[i+4] << 8) | (buf[i+5] << 16);
      if(qq & (1 << 23)){
        // sign-extend
        qq |= 0xff000000;
      }
      std::complex<double> c(std::complex<double>((double) ii, (double) qq));
      one_iq(c);
    }

    if(samples_.size() > rate() * 100){
      fprintf(stderr, "sdrip: samples_ overflow\n");
      samples_.clear();
    }

    samples_mu_.unlock();
  }
}

//
// non-blocking.
// returns all waiting I/Q samples, if any.
// rate_ samples / second.
// sets t to time of last sample.
void
SDRIP::get(std::vector<std::complex<double>> &v)
{
  samples_mu_.lock();
  v.swap(samples_);
  samples_mu_.unlock();
}

//
// read one byte from the TCP socket.
//
int
SDRIP::read8()
{
  char c = 0;

  int cc = read(tcp_, &c, 1);
  if(cc != 1){
    perror("SDRIP TCP socket read");
    exit(1);
  }
  return c & 0xff;
}

//
// read one message from the TCP connection.
//
std::vector<int>
SDRIP::readone()
{
  std::vector<int> v;

  int lsb = read8();
  v.push_back(lsb);
  int msb = read8();
  v.push_back(msb);

  int len = lsb | ((msb & 0x1f) << 8);
  
  if(len < 4 || len > 10000){
    fprintf(stderr, "SDRIP weird len %d, %02x %02x\n", len, lsb, msb);
  }
  //assert(len >= 4 && len < 10000);

  for(int i = 2; i < len; i++){
    int x = read8();
    v.push_back(x);
  }

  return v;
}

std::vector<int>
SDRIP::getitem(int item, int extra)
{
  send_mu_.lock();

  reply_mu_.lock();
  reply_.clear();
  reply_mu_.unlock();
  
  int len = 4;
  if(extra != -1){
    len = 5;
  }
  char buf[5];
  buf[0] = len; // total cmd length
  buf[1] = 0x20; // Request Current Control Item (001 in upper 3 bits)
  buf[2] = item & 0xff; // least-significant byte first
  buf[3] = (item >> 8) & 0xff;
  buf[4] = extra;
  int cc = write(tcp_, buf, len);
  assert(cc == len);

  std::vector<int> v;
  while(1){
    reply_mu_.lock();
    v = reply_;
    reply_mu_.unlock();

    if(v.size() >= 4 && v[2] == (item & 0xff) && v[3] == ((item >> 8) & 0xff)){
      // it's for us!
      break;
    } else {
      usleep(10 * 1000);
    }
  }

  send_mu_.unlock();

  return v;
}

void
SDRIP::setitem(int item, std::vector<int> data)
{
  send_mu_.lock();

  reply_mu_.lock();
  reply_.clear();
  reply_mu_.unlock();

  char *buf = (char*) malloc(data.size() + 4);
  assert(buf);
  
  buf[0] = data.size() + 4; // total cmd length
  buf[1] = 0x00; // Set Current Control Item (000 in upper 3 bits)
  buf[2] = item & 0xff; // least-significant byte first
  buf[3] = (item >> 8) & 0xff;
  for(int i = 0; i < (int) data.size(); i++)
    buf[i+4] = data[i];
  int cc = write(tcp_, buf, data.size()+4);
  assert(cc == (int) data.size()+4);

  std::vector<int> v;
  while(1){
    reply_mu_.lock();
    v = reply_;
    reply_mu_.unlock();

    if(v.size() >= 4 && v[2] == (item & 0xff) && v[3] == ((item >> 8) & 0xff)){
      // it's our's!
      break;
    } else {
      usleep(10 * 1000);
    }
  }

  send_mu_.unlock();

  free(buf);
}
