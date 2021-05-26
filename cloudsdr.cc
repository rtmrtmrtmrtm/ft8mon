//
// RFSpace CloudSDR in cloud mode.
//
// Robert Morris, AB1HL.
//

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
#include "cloudsdr.h"

CloudSDR::CloudSDR()
{
  s_ = -1;
  th_ = 0;
}


void
CloudSDR::open(std::string chan)
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

  s_ = socket(AF_INET, SOCK_STREAM, 0);
  assert(s_ >= 0);

  sockaddr_in sin;
  memset(&sin, 0, sizeof(sin));
  sin.sin_family = AF_INET;
  sin.sin_port = htons(port);
  sin.sin_addr.s_addr = inet_addr(hostname.c_str());

  int ret = connect(s_, (sockaddr *) &sin, sizeof(sin));
  if(ret < 0){
    fprintf(stderr, "CloudSDR: connect(%s) failed\n", hostname.c_str());
    exit(1);
  }

  th_ = new std::thread( [ this ] () { reader(); } );

  // if we send requests too soon, the CloudSDR ignores them.
  usleep(100*1000);

  set_state(1); // stop
  if(hz >= 0){
    set_frequency(hz);
  }
  set_demod(1); // USB
  set_compression(1); // 8 bit 64000 bps raw audio samples
  // set_compression(2); // G711
  set_squelch(-160);
  set_af_filter(100, 3800);
  set_state(2); // run

  printf("name: %s\n", get_name().c_str());
  printf("serial: %s\n", get_serial().c_str());
  printf("status: %s\n", get_status().c_str());
  printf("compression: %d\n", get_compression());
  printf("frequency: %lld\n", get_frequency());
  printf("squelch: %d\n", get_squelch());
  printf("volume: %d\n", get_volume());
  printf("gain: %d\n", get_gain());
  printf("rf filter: %d\n", get_rf_filter());
  {
    int low, high;
    get_af_filter(low, high);
    printf("af filter: %d %d\n", low, high);
  }
  {
    int threshold, slope, decay;
    get_agc(threshold, slope, decay);
    printf("agc: %d %d %d\n", threshold, slope, decay);
  }
}

std::string
CloudSDR::get_name()
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
CloudSDR::get_serial()
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
CloudSDR::get_status()
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
CloudSDR::get_frequency()
{
  std::vector<int> v = getitem(0x0020, 0);

  if((int) v.size() < 10)
    return -1;
  
  long long hz = 0;
  // LSB first
  for(int i = 0; i < 5; i++){
    hz |= ((v[i+5] & 0xff) << (8 * i));
  }
  
  return hz;
}

void
CloudSDR::set_frequency(long long hz)
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

// 0 = LSB
// 1 = USB
// 5 = FM
// 6 = AM
// 8 = WFM
// 10 = DIG
// 11 = RAW
void
CloudSDR::set_demod(int mode)
{
  std::vector<int> data;
  data.push_back(0); // channel ID
  data.push_back(mode);
  setitem(0x0028, data);
}

// 0 = No Output
// 1 = 8 bit 64000 bps raw audio samples
// 2 = 8 bit mu-law G711
void
CloudSDR::set_compression(int mode)
{
  std::vector<int> data;
  data.push_back(0); // channel ID
  data.push_back(mode);
  setitem(0x0084, data);
}

int
CloudSDR::get_compression()
{
  std::vector<int> v = getitem(0x0084, 0);
  if((int) v.size() > 5){
    return v[5];
  } else {
    return -1;
  }
}

// 1 = stop
// 2 = run
void
CloudSDR::set_state(int state)
{
  std::vector<int> data;
  data.push_back(0); // channel ID
  data.push_back(state);
  setitem(0x0018, data);
}

// 0 = forced closed (no audio)
// -160 = forced open (audio always present)
int
CloudSDR::get_squelch()
{
  std::vector<int> v = getitem(0x0080, 0);
  if((int) v.size() < 7)
    return -160;
  int db;
  db = (v[5] & 0xff);
  db |= ((v[6] & 0xff) << 8);
  return db;
}

// 0 = forced closed (no audio)
// -160 = forced open (audio always present)
// seems to default to 0 (no audio)
void
CloudSDR::set_squelch(int sq)
{
  std::vector<int> data;
  data.push_back(0); // channel ID
  data.push_back(sq & 0xff);
  data.push_back((sq >> 8) & 0xff);
  setitem(0x0080, data);
}

int
CloudSDR::get_gain()
{
  std::vector<int> v = getitem(0x0038, 0);
  if((int) v.size() < 6)
    return -1;
  return v[5];
}

// of the headphone output?
// 0..99
int
CloudSDR::get_volume()
{
  std::vector<int> v = getitem(0x0048, 0);
  if((int) v.size() < 6)
    return -1;
  return v[5] & 0xff;
}

int
CloudSDR::get_rf_filter()
{
  std::vector<int> v = getitem(0x0044, 0);
  if((int) v.size() < 6)
    return -1;
  return v[5] & 0xff;
}

// this doesn't match the Rev 0.10 manual:
// there doesn't seem to be an "offset" value.
void
CloudSDR::get_af_filter(int &low, int &high)
{
  std::vector<int> v = getitem(0x0058, 0);
  if((int) v.size() < 7){
    low = high = -1;
    return;
  }
  low = v[5] | (v[6] << 8);
  high = v[7] | (v[8] << 8);
}

void
CloudSDR::set_af_filter(int low, int high)
{
  std::vector<int> data;
  data.push_back(0); // channel ID
  data.push_back(low & 0xff);
  data.push_back((low >> 8) & 0xff);
  data.push_back(high & 0xff);
  data.push_back((high >> 8) & 0xff);
  // but not the "offset" mentioned in the document.
  setitem(0x0058, data);
}

void
CloudSDR::get_agc(int &threshold, int &slope, int &decay)
{
  std::vector<int> v = getitem(0x0050, 0);
  if((int) v.size() < 9){
    threshold = slope = decay = -1;
    return;
  }
  threshold = (signed char) (v[5] & 0xff); // -127 to +128 dB
  slope = v[6]; // 0 to 10 dB
  decay = v[7] | (v[8] << 8); // milliseconds, 20 to 2000
}

void
CloudSDR::reader()
{
  double lastclip = 0;
  
  while(1){
    std::vector<int> v = readone();

    int type = (v[1] >> 5) & 0x7;
    if(type == 0){
      // response to set or request item
      reply_mu_.lock();
      reply_ = v;
      reply_mu_.unlock();
    } else if((int) v.size() == 7 && type == 4){
      // a seven-byte message starting like this is squelched data:
      // 07 80 ...
      // such data starts arriving as soon as a connection is made.
    } else if(type == 4){
      // data
      // 05 82 ss ss type 512 data bytes

      // ssss is signal level, dB * 10.
      // should be -160 to 0 dB (i.e. -1600 to 0).
      // if > 0, overload.
      short ss = (v[2] & 0xff) | ((v[3] & 0xff) << 8);
      if(ss > 0 && now() - lastclip >= 5){
        fprintf(stderr, "CloudSDR CLIP %d\n", ss);
        lastclip = now();
      }

      samples_mu_.lock();
      if(v[4] == 1){
        // not compressed, just 8 bits per sample.
        for(int i = 5; i < (int) v.size(); i++){
          // sign-extend.
          int x = (signed char) (v[i] & 0xff);
          samples_.push_back(x);
        }
      } else if(v[4] == 2){
        // mu-law G711.
        // https://en.wikipedia.org/wiki/G.711
        // each byte: seeeabcd
        //   s is sign
        //   real sample is 1abcd1 left-shifted by eee bits.
        // but some other stuff happens too.
        // yields 14-bit values.
        for(int i = 5; i < (int) v.size(); i++){
          int b = v[i] & 0xff;
          b ^= 0xff;
          int s = ((b & 0x80) != 0); // sign
          int e = ((b >> 4) & 7);    // amount of shift
          int m = (b & 0xf);         // mantissa
          int y = (s ? -1 : 1);
          y *= ((33 + 2*m) * (1 << e) - 33);
          samples_.push_back(y);
        }
      } else {
        fprintf(stderr, "cloudsdr: unknown compression %d\n", v[4]);
      }
      if((int) samples_.size() > 8000 * 100){
        fprintf(stderr, "cloudsdr input overflow\n");
        samples_.clear();
      }
      samples_time_ = now();
      samples_mu_.unlock();
    } else {
      printf("cloudsdr unknown type=%d, ", type);
      for(int i = 0; i < 6 && i < (int) v.size(); i++){
        printf("%02x ", v[i]);
      }
      printf("\n");
    }
  }
}

//
// non-blocking.
// returns all waiting samples, if any.
// always 8000 samples/second.
// sets t to time of last sample.
std::vector<short>
CloudSDR::get(double &t)
{
  std::vector<short> v;
  
  samples_mu_.lock();
  v.swap(samples_);
  t = samples_time_;
  samples_mu_.unlock();

  return v;
}

//
// read one byte from the TCP socket.
//
int
CloudSDR::read8()
{
  char c = 0;

  int cc = read(s_, &c, 1);
  if(cc != 1){
    perror("cloudsdr socket read");
    exit(1);
  }
  return c & 0xff;
}

std::vector<int>
CloudSDR::readone()
{
  std::vector<int> v;

  int lsb = read8();
  v.push_back(lsb);
  int msb = read8();
  v.push_back(msb);

  int len = lsb | ((msb & 0x1f) << 8);
  
  assert(len >= 4 && len < 10000);

  for(int i = 2; i < len; i++){
    int x = read8();
    v.push_back(x);
  }

  return v;
}

std::vector<int>
CloudSDR::getitem(int item, int extra)
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
  int cc = write(s_, buf, len);
  assert(cc == len);

  std::vector<int> v;
  while(1){
    reply_mu_.lock();
    v = reply_;
    reply_mu_.unlock();

    if((int) v.size() >= 4 && v[2] == (item & 0xff) && v[3] == ((item >> 8) & 0xff)){
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
CloudSDR::setitem(int item, std::vector<int> data)
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
  int cc = write(s_, buf, data.size()+4);
  assert(cc == (int) data.size()+4);

  std::vector<int> v;
  while(1){
    reply_mu_.lock();
    v = reply_;
    reply_mu_.unlock();

    if((int) v.size() >= 4 && v[2] == (item & 0xff) && v[3] == ((item >> 8) & 0xff)){
      // it's our's!
      break;
    } else {
      usleep(10 * 1000);
    }
  }

  send_mu_.unlock();

  free(buf);
}

#ifdef CLOUDMAIN
int
main()
{
  CloudSDR sdr;
  sdr.open("192.168.3.140");

  std::vector<double> all;

  while(1){
    sleep(1);
    double ttt;
    std::vector<short> v = sdr.get(ttt);
    printf("%d\n", (int) v.size());
    for(int i = 0; i < (int) v.size(); i++)
      all.push_back(v[i]);
    if((int) all.size() >= 16000){
      writewav(all, "x.wav", 8000);
      writetxt(all, "x.txt");
      printf("wrote x.wav x.txt\n");
      exit(1);
    }
  }
}
#endif
