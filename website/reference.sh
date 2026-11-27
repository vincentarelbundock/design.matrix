#!/usr/bin/env bash
# Convert the Rd files of both packages to Typst reference pages for the
# website: man/*.Rd of xtx into website/reference/, and man/*.Rd of xtx.fit
# into website/reference-fit/ (from the xtx.fit source tree: XTX_FIT_DIR,
# default ../xtx.fit next to this repository). rd2qmd writes plain Typst;
# this adds what a calepin site page needs: the runtime import,
# <website-metadata> with a short sidebar title, and a hidden chunk that
# attaches the package(s) before the examples run.
set -euo pipefail
cd "$(dirname "$0")/.."
fit_dir=${XTX_FIT_DIR:-../xtx.fit}

# convert <man dir> <out dir> <R code of the hidden setup chunk>
convert() {
    local man=$1 out=$2 setup=$3
    rm -rf "$out"
    rd2qmd convert "$man" -o "$out" -f typ -q
    for f in "$out"/*.typ; do
        # sidebar title: the first alias, e.g. "design_matrix"
        alias=$(grep -m1 -o 'aliases: ("[^"]*"' "$f" | sed 's/.*("//; s/"$//')
        {
            printf '#import "/.calepin/calepin.typ" as calepin\n#show: calepin.document\n\n'
            sed -e "s/^#metadata((\$/#metadata((\n  title: \"${alias}\",/" \
                -e "s/^))<rd2qmd>\$/)) <website-metadata>\n\n#calepin.setup(echo: true, eval: true, results: \"verbatim\", warning: false, message: false)\n\n#calepin.chunk(echo: false, results: \"hide\")[\n\`\`\`r\n${setup}\n\`\`\`\n]/" \
                "$f"
        } > "$f.tmp"
        mv "$f.tmp" "$f"
    done
}

convert man website/reference 'library(xtx)'
if [ ! -d "$fit_dir/man" ]; then
    echo "reference.sh: no xtx.fit sources at $fit_dir (set XTX_FIT_DIR)" >&2
    exit 1
fi
convert "$fit_dir/man" website/reference-fit 'library(xtx)\nlibrary(xtx.fit)'
