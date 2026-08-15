CC = gcc
CFLAGS = -fPIC -Wall -O2
LDFLAGS = -shared

all: librscodec.so decoder encoder

librscodec.so: rs_codec.o
	$(CC) $(LDFLAGS) -o $@ $^

rs_codec.o: rs_codec.c rs_codec.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.so

test: test_rs.c librscodec.so
	$(CC) -Wall -O2 -o test_rs test_rs.c -L. -lrscodec -Wl,-rpath,./

decoder: streaming_decoder.c librscodec.so vhs_params.h
	$(CC) -Wall -O2 -o streaming_decoder streaming_decoder.c -lm -lz -L. -lrscodec -Wl,-rpath,./

encoder: encoder.c librscodec.so vhs_params.h
	$(CC) -Wall -O2 -o encoder encoder.c -lm -lz -L. -lrscodec -Wl,-rpath,./
