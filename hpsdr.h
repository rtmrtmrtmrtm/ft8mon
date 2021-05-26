#ifndef hpsdr_h
#define hpsdr_h

#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <math.h>
#include <liquid/liquid.h>
#include <vector>
#include <complex>
#include <mutex>

class HPSDR {
private:
  enum State {
    None,
    Discovering,
    Running,
  };

  State state_;

  // times at which we last did various things.
  double last_discover_;
  double last_general_;
  double last_high_;
  double last_swr_print_;

  // UDP socket.
  int s_;

  int run_;

  // things we learned from the radio.
  sockaddr_in radio_sin_;
  int board_type_;

  // things to tell the radio.
  int do_dither_;
  int do_random_;
  int do_lpf_;
  int do_hpf_;
  int do_attenuate_;
  int do_ext1_; // RX1 from EXT 1 rather than ANT 1
  int do_both_; // for rx, unit 0 on ADC0, units 1-3 on ADC1

  int step_[2]; // attenuation for each of the two ADCs
  double last_clipped_[2]; // time we last increased step_ due to clip
  double last_unclipped_[2]; // time we decreased step_

  // for each of the seven DDCs.
  struct DDC {
    volatile int active_;
    int hz_;
    int sdr_rate_;
    int rate_; // desired rate
    firfilt_crcf filter_;
    long long count_; // of samples, for rate reduction.
    unsigned int expect_seq_;
    std::vector<std::complex<double>> input_; // buffered input, at rate_.
    std::mutex mu_;
  };
  static const int NDDC = 7;
  DDC ddc_[NDDC];

  volatile int configuration_changed_;

  // transmit.
  unsigned int tx_seq_;
  volatile int tx_active_;
  volatile int tx_ptt_;
  volatile int tx_power_;
  volatile int tx_hz_;
  volatile int tx_pa_;
  volatile int tx_drive_;
  volatile int tx_drive_override_;

  static std::vector<char> make_discovery_packet();
  static int open_socket();
  static void discover(int s);
  void tick();
  void check_config_files();
  static int wait_socket(int s, int ms);
  void process();
  void process_discovery_response(const std::vector<char> &buf,
                                  const sockaddr_in &from);
  void send_general();
  void send_ddc();
  void send_duc();
  unsigned int ph(int hz);
  unsigned int make_alex0();
  unsigned int make_alex1();
  unsigned int make_band();
  void send_high();
  void process_packet(const std::vector<char> &buf,
                      const sockaddr_in &from);
  void process_data(const std::vector<char> &buf, int unit);
  void activate(int unit, int rate);
  int watts2drive(int hz, int watts);
  void process_priority_status(const std::vector<char> &buf,
                               const sockaddr_in &from);

  void send_to_duc(const std::vector<std::complex<double>> &iq, int off);

  void loop();
  HPSDR();
  ~HPSDR();
  
 public:
  void lowpass(int unit, double &ii, double &qq);
  
 public:
  static HPSDR *open();
  static void list();
  int allocate_unit(int rate);
  void set_freq(int unit, int hz);
  void get(int unit, std::vector<std::complex<double>> &buf);
  void tx_iq(int hz, const std::vector<std::complex<double>> &v, int rate);
  void tx_real(int hz, const std::vector<double> &v, int rate);
};

#endif
