CXX = clang++ -O
# CXX += -g -fsanitize=address
# CXX = g++9 -O3
FLAGS = -std=c++17 -I/opt/local/include -I/usr/local/include
LIBS = -L/opt/local/lib -L/usr/local/lib -lfftw3 -lsndfile

ft8mon: ft8.cc ft8mon.cc snd.cc libldpc.c osd.cc unpack.cc
	$(CXX) $(FLAGS) ft8mon.cc ft8.cc unpack.cc osd.cc snd.cc libldpc.c -o ft8mon $(LIBS) -lportaudio

clean:
	rm -f ft8mon
