### This file explains how to reproduce the figures of the article
### entitled "Unbiased Hamiltonian Monte Carlo with Couplings", 2017
### by Jeremy Heng and Pierre E. Jacob (Harvard University)
#
# The files of the folder reproduce in debiasedhmc/inst/
# contain the scripts to reproduce the figures.
#
# The files finishing in .run.R are the ones creating the result files;
# these scripts take several hours to run on a desktop computer, even possibly days; you've been warned!
# Note also that the scripts create gigantic files (e.g. ~ 6Gb) so you should
# choose your 'resultsfolder' accordingly!!!!!

# The files finishing in ".plots.R" create the plots, assuming that the files
# finishing in ".run.R" have been run already.
# Indeed the plot files require ".RData" files that are
# created by the ".run.R" files.

# folders
scriptfolder <- "~/path to debiasedhmc/inst/reproduce"
resultsfolder <- "~/path to a folder where large RData files will be created"
#

### Multivariate Normal example (might take a few hours, creates ~ 500 Mb of files)
# run
source(file.path(scriptfolder, "mvnorm.hmctuning.run.R"))
source(file.path(scriptfolder, "mvnorm.contraction.run.R"))
source(file.path(scriptfolder, "mvnorm.meetingtimes.run.R"))
source(file.path(scriptfolder, "mvnorm.meetingtimes.withMH.run.R"))
# plots
source(file.path(scriptfolder, "mvnorm.plots.R"))

### Truncated Multivariate Normal example (takes only a few minutes to run, creates a few Mb of files)
# run
source(file.path(scriptfolder, "tmg.linear.run.R"))
source(file.path(scriptfolder, "tmg.quadratic.run.R"))
# plots
source(file.path(scriptfolder, "tmg.plots"))

### German Credit logistic regression example
# (takes about a day to run, creates ~ 6Gbs of files... !!!)
# run
source(file.path(scriptfolder, "germancredit.hmctuning.run.R"))
source(file.path(scriptfolder, "germancredit.contraction.run.R"))
source(file.path(scriptfolder, "germancredit.meetingtimes.withMH.run.R"))
source(file.path(scriptfolder, "germancredit.meetingtimes.withMH.betterinitrun.R"))
source(file.path(scriptfolder, "germancredit.longtraces.R"))
# plots
source(file.path(scriptfolder, "germancredit.plots.R"))
