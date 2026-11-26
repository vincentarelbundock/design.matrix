#!/usr/bin/env bash
# Convert man/*.Rd to Typst reference pages for the website.
# rd2qmd writes plain Typst; this adds what a calepin site page needs:
# the runtime import, <website-metadata> with a short sidebar title, and a
# hidden chunk that attaches the package before the examples run.
set -euo pipefail
cd "$(dirname "$0")/.."
out=website/reference
rm -rf "$out"
rd2qmd convert man/ -o "$out" -f typ -q
for f in "$out"/*.typ; do
    # sidebar title: the first alias, e.g. "design_matrix"
    alias=$(grep -m1 -o 'aliases: ("[^"]*"' "$f" | sed 's/.*("//; s/"$//')
    {
        printf '#import "/.calepin/calepin.typ" as calepin\n#show: calepin.document\n\n'
        sed -e "s/^#metadata((\$/#metadata((\n  title: \"${alias}\",/" \
            -e 's/^))<rd2qmd>$/)) <website-metadata>\n\n#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)\n\n#calepin.chunk(echo: false, results: "hide")[\n```r\nlibrary(xtx)\n```\n]/' \
            "$f"
    } > "$f.tmp"
    mv "$f.tmp" "$f"
done
