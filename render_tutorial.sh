#!/usr/bin/env bash
# Regenerate vafSOMA_tutorial.html from vafSOMA_tutorial.md (requires pandoc). Usage: bash render_tutorial.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
pandoc vafSOMA_tutorial.md -s --toc --toc-depth=2 --css tutorial_style.css --embed-resources \
  -o vafSOMA_tutorial.html
echo "Wrote vafSOMA_tutorial.html"
