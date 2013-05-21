CC = g++ -O3 -Wall
LIB = ../libim/libim.a ../libjpeg/jpeg/libjpeg.a -lpthread
BIN = jpg2pgm

all: $(BIN)

jpg2pgm: jpg2pgm.cpp
		$(CC) -o jpg2pgm jpg2pgm.cpp $(LIB)

clean:
		/bin/rm -f $(BIN) 
