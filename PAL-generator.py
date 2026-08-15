#!/usr/bin/env python3
import sys
import struct
from PIL import Image

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <sample_rate> <input_image> <output_file.raw>")
    sys.exit(1)

sample_rate = int(sys.argv[1])
image_path = sys.argv[2]
output_file = sys.argv[3]

# Samples per 64 µs line
samples_per_line = int(round(sample_rate * 64e-6))

# Voltage to 8-bit mapping (1.2 V is the max output for my fl2k adapter)
def v_to_u8(v):
    return int(round((v / 1.2) * 255))

SYNC = v_to_u8(0.0)
BLANK_LVL = v_to_u8(0.285)
BLACK_LVL = v_to_u8(0.339)
WHITE = v_to_u8(1.0)

# Active video duration
ACTIVE_DURATION_US = 51.95
active_samples = int(round(sample_rate * ACTIVE_DURATION_US * 1e-6))

def make_samples(duration_us, level):
    #print(f"Duration: {duration_us}, Level: {level}")
    n = int(round(sample_rate * duration_us * 1e-6))
    return [level] * n

def write_line(f, segments, should_pad=False):
    line = []
    for dur, lvl in segments:
        line.extend(make_samples(dur, lvl))
    # pad or trim to exact samples_per_line
    if should_pad:
        if len(line) < samples_per_line:
            line.extend([BLANK_LVL] * (samples_per_line - len(line)))
        elif len(line) > samples_per_line:
            line = line[:samples_per_line]
    f.write(struct.pack(f'{len(line)}B', *line))

def load_image_data(path, target_width, target_height):
    img = Image.open(path).convert('L')  # grayscale
    orig_w, orig_h = img.size

    # 1) Resize keeping aspect ratio so the image fully covers target_height
    scale = max(target_width / orig_w, target_height / orig_h)
    new_w = int(round(orig_w * scale))
    new_h = int(round(orig_h * scale*0.80))
    img = img.resize((new_w, new_h), Image.LANCZOS)

    # 2) Crop to target_height (and target_width) centered
    left = (new_w - target_width) // 2
    top  = (new_h - target_height) // 2
    right = left + target_width
    bottom = top + target_height
    img = img.crop((left, top, right, bottom))

    pixels = img.load()
    data = []
    for y in range(target_height):
        row = []
        for x in range(target_width):
            p = pixels[x, y]  # 0 = black, 255 = white
            # Map to voltage: 0.3V (black) to 1.0V (white)
            v = 0.3 + (p / 255.0) * 0.7
            row.append(v_to_u8(v))
        data.append(row)
    return data

# Total visible lines = 304 * 2 = 608
image_data = load_image_data(image_path, active_samples, 576)

def write_visible_line(f, row_pixels):
    line = []
    # Front porch + sync + back porch
    line.extend(make_samples(4.7, SYNC))
    line.extend(make_samples(5.7, BLANK_LVL))

    # Active video from image
    line.extend(row_pixels)

    # Back porch end
    line.extend(make_samples(1.65, BLANK_LVL))

    # Ensure exact length
    if len(line) < samples_per_line:
        line.extend([BLANK_LVL] * (samples_per_line - len(line)))
    elif len(line) > samples_per_line:
        line = line[:samples_per_line]

    f.write(struct.pack(f'{len(line)}B', *line))

def generate_frame(f, start_line):
    # === FIELD 1 (odd) ===
    # Broad sync (5×)
    for _ in range(5):
        write_line(f, [(27.3, SYNC), (4.7, BLANK_LVL)])

    # Short sync (5×)
    for _ in range(5):
        write_line(f, [(2.35, SYNC), (29.65, BLANK_LVL)])

    for _ in range(17):
        write_line(f, [(4.7, SYNC), (5.7, BLANK_LVL), (51.95, BLACK_LVL), (1.65, BLANK_LVL)], True)

    # Visible lines (odd field)
    for i in range(288):
        write_visible_line(f, image_data[start_line + i * 2 + 1])

    # === FIELD 2 (even) ===
    # Short sync (5×)
    for _ in range(5):
        write_line(f, [(2.35, SYNC), (29.65, BLANK_LVL)])

    # Broad sync (5×)
    for _ in range(5):
        write_line(f, [(27.3, SYNC), (4.7, BLANK_LVL)])

    # Short sync (5×)
    for _ in range(5):
        write_line(f, [(2.35, SYNC), (29.65, BLANK_LVL)])

    # Extra half-line
    write_line(f, [(32, BLANK_LVL)])

    for _ in range(16): #16
        write_line(f, [(4.7, SYNC), (5.7, BLANK_LVL), (51.95, BLACK_LVL), (1.65, BLANK_LVL)], True)

    # Visible lines (even field)
    for i in range(288):
        write_visible_line(f, image_data[start_line + i * 2])

    write_line(f, [(4.7, SYNC), (5.7, BLANK_LVL), (10.0, BLACK_LVL), (10.0, WHITE), (1.6, BLANK_LVL)])

    # Final short sync
    for _ in range(5):
        write_line(f, [(2.35, SYNC), (29.65, BLANK_LVL)])
    print(f.tell(), samples_per_line*625)


with open(output_file, 'wb') as f:
    for frame in range(25):
        generate_frame(f, 0)  # Use top 608 lines; repeat same frame
        print(f"Generated frame {frame + 1}/25")

print(f"Generated 25 frames -> {output_file}")
print(f"samples/line = {samples_per_line}")
print(f"active samples = {active_samples}")
