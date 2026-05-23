#' Compile the `BundibugyoNPI2026.tex` \LaTeX document into a PDF
#'
#' @return Nothing, called for side-effects.
#' @keywords internal
compile_paper <- function() {
  if (!requireNamespace("tinytex", quietly = TRUE)) {
    stop(
      "Package 'tinytex' is required for compile_paper(). ",
      "Install it with install.packages('tinytex').",
      call. = FALSE
    )
  }
  paper_dir <- system.file("paper", package = "BundibugyoNPI2026")
  if (!nzchar(paper_dir)) {
    stop("Could not locate 'paper' directory in installed package.")
  }
  old_wd <- getwd()
  setwd(paper_dir)
  on.exit(setwd(old_wd), add = TRUE)
  tinytex::latexmk("BundibugyoNPI2026.tex")
}
