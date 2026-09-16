#!/usr/bin/env bash
# Regenerate pipeline_subway.svg / .png from pipeline_subway.mmd.
#
# The map is authored in nf-metro's Mermaid-based DSL (see
# https://github.com/seqeralabs/nf-metro) and rendered with the real
# nf-metro CLI, so it matches the same nf-core subway/metro-map style used
# for pipelines like nf-core/rnaseq, rather than a hand-drawn approximation.
#
# This cluster's system Python is externally managed (PEP 668), so nf-metro
# is installed into a throwaway virtualenv rather than system-wide/--user.
# PNG export goes through cairosvg since no system rsvg-convert/inkscape is
# available here; swap in `rsvg-convert` or `inkscape` if you have one.
set -euo pipefail
cd "$(dirname "$0")"

VENV=.render-venv
if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
    "$VENV/bin/pip" install --quiet "nf-metro==1.1.0" cairosvg
fi

# --no-chrome-css: bake colors as static presentation attributes instead of
# CSS var()s, which cairosvg (and some other rasterizers/PowerPoint's SVG
# import) cannot resolve.
"$VENV/bin/nf-metro" render pipeline_subway.mmd \
    -o pipeline_subway.svg \
    --theme nfcore-light \
    --embed-font \
    --no-chrome-css

"$VENV/bin/python" -c "
import cairosvg
cairosvg.svg2png(url='pipeline_subway.svg', write_to='pipeline_subway.png', output_width=3000)
"

echo "Wrote pipeline_subway.svg and pipeline_subway.png"
echo "(delete $VENV/ any time; it is regenerated on the next run)"
