# dartR
Importing and Analysing DArT type snp and silicodart data


Currently the installation of the package might not run smoothly as it requires additional bioconductor and github package. 
A that always seems to work is copy paste the script below into your R-console. If you are lucky it simply installs all packages:


install.packages("devtools")
library(devtools)

source("http://bioconductor.org/biocLite.R")
biocLite("qvalue", suppressUpdates=T)
biocLite("SNPRelate", suppressUpdates=T)
install_github("whitlock/OutFLANK")
install_github("green-striped-gecko/dartR")

library(dartR)



Unfortunately sometime to us unkown reason R is not able to instlal app dependent packages and breaks with an error message. 
Then you need to install packages "by hand". For example you may find:

ERROR: dependency 'seqinr' is not available for package 'dartR'
* removing 'C:/Program Files/R/library/dartR'
Error: Command failed (1)

Then you need to install the package seqinr via: 
install.packages("seqinr")

And run the last line of code again:
install_github("green-striped-gecko/dartR")

This "game"  of install.package() and install_github() [the last two steps] might need to be repeated for additional packages
as for whatever reason R once broken does no longer install all packages (if anyone could tell me a way to fix this it would be highly 
appreciated). Finally you should be able to run:

install_github("green-striped-gecko/dartR")
library(dartR)

without any error and you are done. 

Have fun working with dartR (any issues your encounter please use the issues tab on github or via email to me or Arthur via 
glbugs@aerg.canberra.edu.au

Cheers, Bernd & Arthur


