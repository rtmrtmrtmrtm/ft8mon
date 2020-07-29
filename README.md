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

The columns are HHMMSS, SNR, time offset, frequency, and the message.

To read input from a recorded WAV file:

```
  ./ft8mon -file xxx.wav
```

If you have an Airspy HF+ Discovery, and you install the airspyhf and
liquid dsp libraries, you can uncomment the relevant lines in the
Makefile.

Robert Morris,
AB1HL
