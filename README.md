# ft8mon
Demodulate the FT8 WSJT-X protocol of Taylor and Franke.
Input from a sound card via portaudio.
Written in C++, using FFTW.
Runs on Linux, MacOS, and FreeBSD.

To compile:

```
  make
```

To get a list of sound card numbers:

```
  ./ft8mon -list
```

To listen to sound card X, left channel:

```
  ./ft8mon -card X 0
```

You should see output like this:

```
094445  -6 -0.7 1083 JA2KGH CM5VVC EL92
094445 -15  0.1 2584 FG5GH   N4EFS R-07
094445 -10  0.3  601 WA4PT  VK5MRD RRR 
094445 -19  0.0  820 CQ   DX WD6DBM CM97
094445 -15 -0.0 2398 CQ VK3HGQ QF33
094445 -15 -0.4 1227 AB1HL  NA7K   -10
```

The columns are HHMMSS, SNR, time offset, audio frequency, and the message.

To read input from a recorded WAV file:

```
  ./ft8mon -file xxx.wav
```

For Airspy HF+ Discovery support, install the airspyhf
and liquid dsp libraries, and uncomment the relevant lines in the
Makefile. For RFspace SDR-IP, NetSDR, CloudIQ, and CloudSDR
support, install liquid dsp and edit the Makefile. Similarly
for the Apache ANAN-7000dle. Then try commands like these:

```
  ./ft8mon -card airspy ,14.074
  ./ft8mon -card hpsdr 192.168.3.100,14.074
  ./ft8mon -card sdrip 192.168.3.100,14.074
```

Robert Morris,
AB1HL
