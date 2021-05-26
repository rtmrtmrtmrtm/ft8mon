#ifndef sdrip_h
#define sdrip_h

#include <thread>
#include <mutex>
#include <liquid/liquid.h>
#include <vector>
#include <complex>
#include <string>

class SDRIP {
private:
  int rate_; // desired output rate
  int iq_rate_; // SDRIP delivers at this rate
  int tcp_; // TCP socket, for control.
  int udp_; // UDP socket, to receive I/Q data;
  int adc_rate_; // to distinguish SDRIP/NetSDR from Cloud*

  // TCP reader thread
  std::thread *control_th_;

  // only one writer (for e.g. set_frequency()) at a time.
  std::mutex send_mu_;

  // control_reader puts latest reply here.
  std::vector<int> reply_;
  std::mutex reply_mu_; // reply_

  // UDP reader thread
  std::thread *data_th_;

  // written by data_reader, read by get().
  // I/Q at rate_.
  std::vector<std::complex<double>> samples_;
  std::mutex samples_mu_;

  // iq_rate_ -> rate_ conversion
  firfilt_crcf filter_;
  long long iqcount_;
  long long outcount_;

  int read8();
  std::vector<int> readone();
  void control_reader();
  void data_reader();
  std::vector<int> getitem(int item, int extra = -1);
  void setitem(int item, std::vector<int> data);
  void one_iq(const std::complex<double> iq);

public:
  SDRIP(int rate);
  void open(std::string);
  std::string get_name();
  std::string get_serial();
  std::string get_status();
  long long get_frequency();
  void set_frequency(long long hz);
  void set_state(int state);
  void set_iq_rate(int rate);
  int get_iq_rate();
  int get_gain();
  void set_gain(int gain);
  int get_ad_mode();
  void set_ad_mode(int mode);
  int get_rf_filter();
  void get_udp_addr(unsigned int &ip, int &port);
  void set_udp_addr(unsigned int ip, int port);
  int get_adc_rate();
  int get_packet_size();
  void set_packet_size(int sz);
  bool is_sdrip() { return adc_rate_ > 70000000 && adc_rate_ < 90000000; }
  bool is_cloudsdr() { return adc_rate_ > 120000000 && adc_rate_ < 130000000; }
  void get(std::vector<std::complex<double>> &);
  int rate() { return rate_; }
};

#endif
