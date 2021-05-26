CXX = c++ -O
# CXX += -g -fsanitize=address
# CXX = g++9 -O3
FLAGS = -std=c++17 -I/opt/local/include -I/usr/local/include -I/opt/local/include/libairspyhf -I/usr/include/libairspyhf
LIBS = -L/opt/local/lib -L/usr/local/lib -lfftw3 -lsndfile

MOREC = 
MOREH = 

# uncomment if you have the airspyhf and liquid dsp libraries.
# try -lusb or -lusb-1.0 ; also apt install libusb-1.0-0-dev
# CXX += -DUSE_AIRSPYHF
# LIBS += -lairspyhf -lliquid -lusb-1.0

# for the Apache ANAN-7000dle, and possibly other HPSDR radios.
# CXX += -DUSE_HPSDR
# MOREC += hpsdr.cc
# MOREH += hpsdr.h
# LIBS += -lliquid

# for the RFSpace SDR-IP, NetSDR, CloudIQ and CloudSDR in I/Q mode.
# CXX += -DUSE_SDRIP
# MOREC += sdrip.cc
# MOREH += sdrip.h
# LIBS += -lliquid

ft8mon: ft8.cc ft8mon.cc snd.cc libldpc.c osd.cc unpack.cc util.cc fft.cc cloudsdr.h cloudsdr.cc $(MOREC) $(MOREH)
	$(CXX) $(FLAGS) ft8mon.cc ft8.cc unpack.cc osd.cc snd.cc util.cc fft.cc libldpc.c cloudsdr.cc $(MOREC) -o ft8mon $(LIBS) -lportaudio -pthread

clean:
	rm -f ft8mon
