options(repos = unique( c(
  mvb = "https://markbravington.r-universe.dev",
  getOption( "repos")[ "CRAN"])))
install.packages( "debug")