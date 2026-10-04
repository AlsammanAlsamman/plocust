#!/bin/bash
# Render flowchart.html to ../images/flowchart.png with headless Chrome at 2x, then crop to the 1800x960 chart.
cd "$(dirname "$0")"
OUT="$PWD/../images/flowchart.png"
google-chrome --headless=new --disable-gpu --hide-scrollbars --allow-file-access-from-files \
  --force-device-scale-factor=2 --window-size=1800,1200 --virtual-time-budget=3000 \
  --screenshot="$OUT" "file://$PWD/flowchart.html" 2>/dev/null
python3 - "$OUT" <<'PY'
import sys
from PIL import Image
im = Image.open(sys.argv[1]); im.crop((0, 0, 3600, 1920)).save(sys.argv[1], optimize=True)
PY
ls -la "$OUT"
