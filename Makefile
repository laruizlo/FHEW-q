PREFIX=$(HOME)
CC=g++
AR=ar
CFLAGS= -O3 -Wall -fopenmp
LDFLAGS= -L. -lfhemb -lfftw3 -lfftw3l
INCLUDE=distrib.h tools.h gates.h LWE.h HomACC.h FFT.h params.h

all: libfhemb.a genFHEW

#cmd: cmd/gen cmd/enc cmd/nand cmd/dec

install: $(INCLUDE) libfhemb.a
	install $(INCLUDE) $(PREFIX)/include
	install libfhemb.a $(PREFIX)/lib

uninstall:
	rm $(PREFIX)/lib/libfhemb.a
	rm $(PREFIX)/lib/{distrib,tools,gates,LWE,HomACC,FFTRing,params}.h

clean:
	rm *.o libfhemb.a genFHEW || echo nothing to clean

libfhemb.a: distrib.o tools.o gates.o FFTRing.o LWE.o HomACC.o time_profiler.o FHEtools.o
	$(AR) -q libfhemb.a distrib.o tools.o gates.o FFTRing.o LWE.o HomACC.o time_profiler.o FHEtools.o

distrib.o: distrib.cpp distrib.h params.h
	$(CC) $(CFLAGS) -c distrib.cpp

tools.o: tools.h tools.cpp
	$(CC) $(CFLAGS) -c tools.cpp

gates.o: gates.h gates.cpp params.h LWE.h HomACC.h
	$(CC) $(CFLAGS) -c gates.cpp

FFTRing.o: FFTRing.h FFTRing.cpp params.h HomACC.h
	$(CC) $(CFLAGS) -c FFTRing.cpp

LWE.o: LWE.h LWE.cpp FFTRing.h params.h distrib.h
	$(CC) $(CFLAGS) -c LWE.cpp

HomACC.o: HomACC.h HomACC.cpp FFTRing.h LWE.h params.h
	$(CC) $(CFLAGS) -c HomACC.cpp

time_profiler.o: time_profiler.h time_profiler.cpp
	$(CC) $(CFLAGS) -c time_profiler.cpp

genFHEW: genFHEW.cpp libfhemb.a
	$(CC) $(CFLAGS) -o genFHEW genFHEW.cpp $(LDFLAGS)



#common.o: cmd/common.cpp cmd/common.h libfhemb.a
#	$(CC) $(CFLAGS) -c cmd/common.cpp 

#cmd/gen: cmd/gen.cpp common.o libfhemb.a
#	$(CC) $(CFLAGS) -o cmd/gen cmd/gen.cpp common.o $(LDFLAGS)

#cmd/enc: cmd/enc.cpp common.o libfhew.a
#	$(CC) $(CFLAGS) -o cmd/enc cmd/enc.cpp common.o $(LDFLAGS)

#cmd/nand: cmd/nand.cpp common.o libfhew.a
#	$(CC) $(CFLAGS) -o cmd/nand cmd/nand.cpp common.o $(LDFLAGS)

#cmd/dec: cmd/dec.cpp common.o libfhew.a
#	$(CC) $(CFLAGS) -o cmd/dec cmd/dec.cpp common.o $(LDFLAGS)
