#!/bin/bash
# Render atlas.html to ../images/rice_atlas.png (headless Chrome at 2x, light theme), cropped to the 2400x1400 figure.
cd "$(dirname "$0")"
OUT="$PWD/../images/rice_atlas.png"
google-chrome --headless=new --disable-gpu --hide-scrollbars --allow-file-access-from-files \
  --force-device-scale-factor=2 --window-size=2400,1800 --virtual-time-budget=4000 \
  --blink-settings=preferredColorScheme=1 \
  --screenshot="$OUT" "file://$PWD/atlas.html" 2>/dev/null
python3 - "$OUT" <<'PY'
import sys
from PIL import Image
im = Image.open(sys.argv[1]); im.crop((0, 0, 4800, 3040)).save(sys.argv[1], optimize=True)
PY
ls -la "$OUT"
