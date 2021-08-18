# set up renv for reproducibility

if (!requireNamespace("remotes"))
  install.packages("remotes")

remotes::install_github("rstudio/renv")


# initialize project-specific saved environment
renv::init()


# whenever packages installed or removed:
renv::snapshot()


# when transferring project to new machine, start with this to install project specific environment
renv::restore()
