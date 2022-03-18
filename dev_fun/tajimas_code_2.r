### Set my working directory
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the sample data file
my.data = read.vcf("heli_phased.vcf", header=TRUE, stringsAsFactors=FALSE)

### Divide the data into windows, and set up
### a results table
my.win=seq(min(my.data$POS), max(my.data$POS), by=600) # This command gets all of the window Start positions
my.results=data.frame(Start=my.win, End=(my.win+600), SNPs=rep(0, length(my.win))) # By adding 600 to every Start, I get all of the Window End positions

### Calculate the number of SNPs in each window
for (i in 1:nrow(my.results)) { # Notice here that the number of iterations in my loop is equal to the number of rows in my RESULTS, not my input data
  d=subset(my.data, (my.data$POS>=my.results$Start[i] & my.data$POS<my.results$End[i])) # For each of my windows, I get the subset of data that falls into that window
  my.results$SNPs[i]=nrow(d) # Then I count the number of SNPs in that subset and save it in my results table
}

### Fill in the Number of Chromosomes for each Window
num.samples=ncol(my.data)-9
my.results$Nchr=2*(num.samples) # Because there is no missing data in this particular file, I can assume the value for 2N is the same for every window

### Calculate Waterson's Theta for Each window
a=seq(from=1, to=((2*num.samples)-1), by=1) # Like our Theta calculations last week, I need to get the sequence of numbers from 1 to 2N-1,
my.results$ThetaW = my.results$SNPs/(sum(1/a)) # then theta-w = #SNPs/(1/1 + 1/2 + 1/3 + ... + 1/(2N-1))

### Create a function for derived counts
derivedCount <- function(row) {
  row=as.vector(row, mode="character")
  x=count.genotypes(get.field(row[10:length(row)], row[9], "GT"))
  dc=(2*x["aa"])+x["Aa"]
  return(unname(dc))
}

### Now, calculate Pi for each window
my.results$Pi=rep(0, nrow(my.results)) # Here, I add a column of zeros to my results table so I can fill in the Pi calculations
for (i in 1:nrow(my.results)) { # Again, I want a calculation for every WINDOW, not every SNP in the VCF file, so my loop goes through 1:50 (there are 50 windows)
  d=subset(my.data, (my.data$POS>=my.results$Start[i] & my.data$POS<my.results$End[i])) # get the subset in the window
  j=apply(d, 1, FUN=derivedCount) # The apply function gets the derived count at EVERY row in my subset of data; so the "j" variable is a VECTOR of counts
  # this is the same as doing a second loop inside of the first loop; it is just a little more efficient
  c=rep(my.results$Nchr[i], length(j)) # I don't actually have to do this, but here I am make a vector of my 2N values (same as c in Pi equation) that is the same length as my vector of j values
  my.results$Pi[i]=sum((2*j*(c-j))/(c*(c-1))) # Finally, I use the equation for Pi to get the calculation for the whole window
}

### Calculate Var(d) for each window
my.results$varD=rep(0, nrow(my.results)) # Set up an empty column to hold the results of the variance function

### Create a function to calculate var(d)
### n = the number of chromosomes (2 x num. samples)
### S = # SNPs
variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}

for (i in 1:nrow(my.results)) { # Loop through the windows one more time, this time calculate the variance for each window
  my.results$varD[i] = variance.d(n=my.results$Nchr[i], S=my.results$SNPs[i])
}

### Now, calculate Tajima's D
my.results$TajimaD = (my.results$Pi - my.results$ThetaW)/(sqrt(my.results$varD)) # With all of the calculations set up in my table, I can do the last calculation for Tajima's D using R's vector math capabilities (and trust that it understands I want one calculation per row in the table)
plot(my.results$Start, my.results$TajimaD, pch=20,xlab="Position", ylab="Tajima's D", type="l", ylim=c(-2.2,0)) 

### Add rectangles at the top to show genes
rect(68023734,-0.4,68024960,-0.3, col="goldenrod")
rect(68028443,-0.4,68029634,-0.3, col="goldenrod")
rect(68034103,-0.4,68036074,-0.3, col="goldenrod")

### Add a rectangle for positive selection zone
rect(68010000, -3, 68060000, -2, col=rgb(1,0,0,0.4), border=NA)

