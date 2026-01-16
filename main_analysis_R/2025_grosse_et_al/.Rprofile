source("renv/activate.R")

message("renv status")
renv::status()
message()

## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}

# disable omnipathr writing logfiles to the repo basepath
# did not find any documentation on how to change the omnipathr logging path itself
options(omnipathr.logfile='none')

