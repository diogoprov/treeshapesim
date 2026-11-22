# autotest traces the package with typetracer, whose code injected into each
# function draws random numbers, so trees frozen for a seed cannot be matched
skip_if_traced <- function() {
  testthat::skip_if(nzchar(Sys.getenv("TYPETRACER_LEAVE_TRACES")),
                    "typetracer draws random numbers while tracing")
}
