CC = g++ -O3 -Wall

all: libim.a

libim.a: im*.cpp im*.h
	$(CC) -c im_color.cpp im_complex.cpp im_float.cpp im_short.cpp im_pipeline.cpp im_queue.cpp
	ar rc libim.a *.o
	/bin/rm -f *.o

clean:
	/bin/rm -f *.o libim.a