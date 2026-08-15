# Store data on VHS
Store digita data on VHS tapes using software and fl2k (DAC) and cxadc (ADC)

## Demo
Here is a [Video](https://mango.vg/files/vhs.mp4) of RAW PCM audio (CD quality) being played directly from tape. \
At this speed on my setup the error correction is not perfect, but good enough for PCM. \
80-100 bits per line is way more resilient.

## Build
```bash
git clone https://github.com/IlyaMZP/vhs_data_store.git
cd vhs_data_store
make
```

## Run
```bash
# Encode
./encoder ./encoder.c output.u8
# Decode
cat output.u8 | ./streaming_decoder - output.txt
```

Currently the encoder doesn't pad the useful data, so direct round trip will fail. \
You can generate a pad file with `PAL-generator.py` with some image and simply concatenate them. \
Because `streaming_decoder` only calibrates levels in the beginning, you should use an image with good contrast.

## Writing tapes
```bash
mkfifo write_pipe
./encoder ./encoder.c write_pipe
cat pad.u8{,,,} write_pipe | ncat -lkp 1234
sudo ~/osmo-fl2k/build/src/fl2k_tcp -s 40000000 # You'll need to patch fl2k_tcp to support uint8_t samples
```

## Reading tapes
```bash
# Assuming your CX card has 40MHz crystal
pv /dev/cxadc0 | ./streaming_decoder - out.bin
```

You might have to play with `vhs_params.h`

## Word of caution
I didn't have time to finish this project, but really wanted to try and store data on VHS tapes. The implementation is quite awful, since I used a clanker to combine my snippets into a finished piece. \
I do not condone the use of AI, but without it this silly experiment would not have been finished.
