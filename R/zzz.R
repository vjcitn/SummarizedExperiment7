# Required by S7 so that methods registered for base/internal generics
# (dim, dimnames, length, names, `[`, rbind, cbind, print) are activated
# when the package is loaded normally via library(), not just when its
# source is run directly in the current session (e.g. devtools::load_all()).
# See ?S7::methods_register.

.onLoad <- function(...) {
  S7::methods_register()
}
