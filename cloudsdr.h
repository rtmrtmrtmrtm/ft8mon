#ifndef cloudsdr_h
#define cloudsdr_h

#include <thread>
#include <mutex>
#include <vector>
#include <string>

class CloudSDR {
private:
  int s_;

  // reader thread
  std::thread *th_;

  // only one writer (for e.g. set_frequency()) at a time.
  std::mutex send_mu_;

  // reader puts latest reply here.
  std::vector<int> reply_;
  std::mutex reply_mu_; // reply_

  // written by reader, read by get().
  std::vector<short> samples_;
  double samples_time_; // when last sample arrived.
  std::mutex samples_mu_;

  int read8();
  std::vector<int> readone();
  void reader();
  std::vector<int> getitem(int item, int extra = -1);
  void setitem(int item, std::vector<int> data);

public:
  CloudSDR();
  void open(std::string);
  std::string get_name();
  std::string get_serial();
  std::string get_status();
  long long get_frequency();
  void set_frequency(long long hz);
  void set_demod(int mode);
  int get_compression();
  void set_compression(int mode);
  void set_state(int state);
  int get_squelch();
  void set_squelch(int sq);
  int get_volume();
  int get_gain();
  int get_rf_filter();
  void get_af_filter(int &low, int &high);
  void set_af_filter(int low, int high);
  void get_agc(int &threshold, int &slope, int &decay);
  std::vector<short> get(double &);
  int rate() { return 8000; }
};

#endif
