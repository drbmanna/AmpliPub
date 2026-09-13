# Renders report.Rmd for the workflow's `report` rule. Knitting happens in a temporary
# directory, so nothing is written into the repository checkout.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

out <- normalizePath(snakemake@output[[1]], mustWork = FALSE)
work <- tempfile("amplipub-report-")
dir.create(work)

rmarkdown::render(
  input = file.path(snakemake@scriptdir, "report.Rmd"),
  output_file = basename(out),
  output_dir = dirname(out),
  intermediates_dir = work,
  knit_root_dir = work,
  params = list(
    results = normalizePath(snakemake@input[["results"]]),
    run_info = normalizePath(snakemake@input[["run_info"]]),
    provenance_dir = normalizePath(dirname(snakemake@input[["environments"]]))
  ),
  envir = new.env(),
  quiet = FALSE
)
