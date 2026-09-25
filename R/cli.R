## Helpers for the command-line scripts in inst/scripts.

## Parse "--key value" and "--key=value" arguments over a list of defaults.
## Dashes in keys become underscores, so --max-width sets max_width.  A bare
## flag sets its key to "true".  Unknown keys are an error, which catches typos
## before a long job starts.
fc_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  a <- defaults; i <- 1L
  while (i <= length(args)) {
    x <- args[i]
    if (!startsWith(x, "--")) stop("unexpected argument '", x, "'")
    x <- sub("^--", "", x)
    if (grepl("=", x, fixed = TRUE)) {
      k <- sub("=.*$", "", x); v <- sub("^[^=]*=", "", x); i <- i + 1L
    } else if (i < length(args) && !startsWith(args[i + 1L], "--")) {
      k <- x; v <- args[i + 1L]; i <- i + 2L
    } else {
      k <- x; v <- "true"; i <- i + 1L
    }
    k <- gsub("-", "_", k)
    if (!k %in% names(defaults))
      stop("unknown option --", k, "; known: ",
           paste0("--", gsub("_", "-", names(defaults)), collapse = " "))
    a[[k]] <- v
  }
  a
}

## Stop with a usage line if a required option is empty.
fc_require <- function(a, keys, usage) {
  miss <- keys[!nzchar(unlist(a[keys]))]
  if (length(miss)) stop("missing ", paste0("--", gsub("_", "-", miss), collapse = ", "),
                         "\nusage: ", usage, call. = FALSE)
  invisible(TRUE)
}

## Path to a script shipped with the package.
fc_script <- function(name = NULL) {
  d <- system.file("scripts", package = "fragcorr", mustWork = TRUE)
  if (is.null(name)) d else file.path(d, name)
}

## Posterior files in a directory, read into a named list (names from the file
## names), as region_tables() and decay_curves() take them.
read_posteriors <- function(dir) {
  files <- list.files(dir, "[.]rds$", full.names = TRUE)
  if (!length(files)) stop("no posterior .rds files in ", dir)
  out <- lapply(files, function(f) readRDS(f)$posterior)
  names(out) <- sub("[.]rds$", "", basename(files))
  out
}
