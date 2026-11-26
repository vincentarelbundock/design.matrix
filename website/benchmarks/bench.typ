// Show benchmark results saved by `make bench` (website/benchmarks/run.R).
// Building the site never runs a benchmark: it reads results/<name>.csv.

#import "/.calepin/calepin.typ" as calepin

// Display an R script without running it: calepin executes every raw block
// whose language is an engine, so turn that off around this one.
#let bench-code(name) = {
  calepin.setup(fenced-chunks: false)
  raw(read("/benchmarks/" + name + ".R"), lang: "r", block: true)
  calepin.setup(fenced-chunks: true)
}

#let bench-table(name, code-column: true) = {
  let rows = csv("/benchmarks/results/" + name + ".csv")
  let ncol = rows.at(0).len()
  let first(c) = if code-column { raw(c) } else { [#c] }
  table(
    columns: ncol,
    align: (left,) + (right,) * (ncol - 1),
    table.header(..rows.at(0).map(h => strong(h))),
    ..rows.slice(1).map(r => (first(r.at(0)),) + r.slice(1).map(c => [#c.trim()])).flatten(),
  )
}

// The benchmark's code, then its results.
#let bench(name) = {
  bench-code(name)
  bench-table(name)
}
