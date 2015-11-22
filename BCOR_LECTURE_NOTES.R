## ----tidy=TRUE-----------------------------------------------------------
rm(list=ls())
p0 <- 0.95
t <- 1:100000

marks <- c(1,50000,100000)
u <- 0.000001
qt <- 1 - p0*exp(-u*t)
qt5 <- 1 - p0*exp(-0.00001*t)
qt4 <- 1 - p0*exp(-0.0001*t)

plot(x=t, y=qt, xlab="Time (Generations)", ylab="q (frequency of B allele)", type="l", col="red", xaxt="n", ylim=c(0,1))
axis(1,at=marks,labels=marks)
points(x=t,y=qt5,type="l",col="blue")
points(x=t,y=qt4,type="l",col="orange")

## ----tidy=TRUE, echo=FALSE-----------------------------------------------

# FUNCTION for plotting results of common garden experiments
# 21 November 2015
# NJG

CommonGarden <- function(TraitMeans=c(40,20,70)){
  Pop1 <- rnorm(15,TraitMeans[1],10)
  Pop2 <- rnorm(15,TraitMeans[2],10)
  Pop3 <- rnorm(15,TraitMeans[3],10)
PopData <- c(Pop1,Pop2,Pop3)
Treatment <- rep(c("P1","P2","P3"),each=15)
par(mar=c(6,6,4,2))
BoxPlotDataFrame <- data.frame(x=Treatment,y=PopData)
boxplot(y~x,data=BoxPlotDataFrame,ylab="Trait Value",sub="Donor Source",cex.sub=2,cex.axis=1.5, cex.lab=2,col="bisque")
  }

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
CommonGarden(c(40,20,70))

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
CommonGarden(c(35,35,35))


## ----tidy=TRUE, echo=FALSE-----------------------------------------------
RecipPlot <- function(ydata=c(10,20,30,5)){
  # set margins and set up axis locations
  par(mar=c(9,4,4,2))
  xdata <- c(1,1,2,2)
  
  # set up an empty plot
  plot(x=xdata,y=ydata,xlim=c(0.5,2.5),ylim=c(min(ydata)-5,max(ydata)+5),ann=F,axes=F,type="n")
  grid(nx=0,ny=10)
  # add x axis and labels
  axis(side=1,labels=c("Cold","Warm"),at=c(1,2),tick=T, cex.axis=1.5)
  mtext("Donor Population",side=1,cex=2,line=4)
  mtext("Trait Value",side=2,cex=2, line=1)
  box()
  
  # add lines and points
  lines(c(1,2),ydata[c(1,2)])
  lines(c(1,2),ydata[c(3,4)])
  points(c(1,2),ydata[c(1,2)],cex=4,pch=21,bg="skyblue")
  points(c(1,2),ydata[c(3,4)],cex=3,pch=22,bg="salmon")
  
  # add legend
  legend("topright",legend=c("Cold Recipient Site","Warm Recipient Site"),pch=c(21,22),pt.cex=2,pt.bg=c("skyblue","salmon"))
  
}
RecipPlot(c(10,10,10,10))

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
RecipPlot(c(20,20,10,10))

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
RecipPlot(c(20,10,20,10))

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
RecipPlot(c(20,5,25,10))

## ----tidy=TRUE, echo=FALSE-----------------------------------------------
RecipPlot(c(20,50,15,10))

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate observed allele frequencies or a single gene with 2 alleles  
# NJG
# 21 November 2015

AlleleFreq_2A <- function(x=c(AA=100, AB=50, BB=50)) {
  
  # Pull out counts of individual genotypes
  AA <- x[1]
  AB <- x[2]
  BB <- x[3]
  # Create a vector and divide by the sum; works for frequencies or raw counts as input
  Gen_Freq <- c(AA,AB,BB)/sum(AA + AB + BB)
  
  # Print genotype frequencies
  cat("Observed genotypic frequencies:", "\n","freq(AA) = ", Gen_Freq[1], "\n", "freq(AB) = ", Gen_Freq[2], "\n", "freq(BB) = ", Gen_Freq[3],"\n")
  cat("\n")
  
  # Create vector for allele frequencies
  Allele_Freq <- vector("numeric",2)
  
  # Use genotypes to calculate allele frequencies
  Allele_Freq[1] <- Gen_Freq[1] + 0.5*Gen_Freq[2]
  Allele_Freq[2] <- Gen_Freq[3] + 0.5*Gen_Freq[2]
  
  # Print allelic frequencies
  cat("Observed allelic frequencies:", "\n", "freq(A) = ", Allele_Freq[1], "\n", "freq(B) = ", Allele_Freq[2], "\n")
  cat("\n")
  # Return the output vector
  return(Allele_Freq)
}

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate observed allele frequencies or a single gene with 3 alleles  
# NJG
# 21 November 2015

AlleleFreq_3A <- function(x=c(JJ=100, JK=50, JL=50, KL=50, KK=50, LL=50)) {
  
  # Convert input vector to individual genotypes
  JJ <- x[1]
  JK <- x[2]
  JL <- x[3]
  KL <- x[4]
  KK <- x[5]
  LL <- x[6]
  # Create a vector and divide by the sum; works for frequencies or raw counts as input
  Gen_Freq <- c(JJ, JK, JL, KL, KK, LL)/sum(JJ, JK, JL, KL, KK, LL)
  
  # Print genotype frequencies
  cat("Observed genotypic frequencies:", "\n","freq(JJ) = ", Gen_Freq[1], "\n", "freq(JK) = ", Gen_Freq[2], "\n", "freq(JL) = ", Gen_Freq[3],"\n", "freq(KL) = ", Gen_Freq[4], "\n", "freq(KK) = ", Gen_Freq[5], "\n", "freq(LL) = ", Gen_Freq[6],"\n")
  cat("\n")
  
  # Create vector for allele frequencies
  Allele_Freq <- vector("numeric",3)
  
  # Use genotypes to calculate allele frequencies
  Allele_Freq[1] <- Gen_Freq[1] + 0.5*Gen_Freq[2] + 0.5*Gen_Freq[3]
  Allele_Freq[2] <- Gen_Freq[5] + 0.5*Gen_Freq[2] + 0.5*Gen_Freq[4]
  Allele_Freq[3] <- Gen_Freq[6] + 0.5*Gen_Freq[3] + 0.5*Gen_Freq[4]
  
  # Print allelic frequencies
 cat("Observed allelic frequencies:", "\n", "freq(J) = ", Allele_Freq[1], "\n", "freq(K) = ", Allele_Freq[2], "\n", "freq(L) = ", Allele_Freq[3], "\n")
  cat("\n")
 
  # Return the output vector
  return(Allele_Freq)
}

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate Hardy-Weinberg genotypic frequency for a single gene with 2 alleles  
# NJG
# 21 November 2015

HardyWeinberg_2A <- function(x=c(p=0.7, q=0.3)){
  
  # Convert input vector to individual frequencies
  p <- x[1]
  q <- x[2]
# Create a vector for genotypic frequencies
Genotype_Freq <-vector("numeric",3)  

# Use Hardy-Weinberg equation to calculate genotypic frequencies from allelic frequencies
Genotype_Freq[1] <- p^2
Genotype_Freq[2] <- 2*p*q
Genotype_Freq[3] <- q^2

# Print allelic frequencies
cat("Observed allelic frequencies:", "\n", "f(A) = ", p, "\n", "f(B) = ", q, "\n")
cat("\n")

# Print expected Hardy-Weinberg genotypic frequencies
cat("Expected Hardy-Weinberg genotypic frequencies:", "\n", "H-W f(AA) = ", Genotype_Freq[1], "\n", "H-W f(AB) = ", Genotype_Freq[2], "\n", "H-W f(BB) = ", Genotype_Freq[3], "\n")
  cat("\n")

  # Return the output vector
  return(Genotype_Freq)
}

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate Hardy-Weinberg genotypic frequency for a single gene with 3 alleles  
# NJG
# 21 November 2015

HardyWeinberg_3A <- function(x=c(p=0.7, q=0.2, r=0.1)){
  
  # Convert input vector into individual allelic frequencies
  p <- x[1]
  q <- x[2]
  r <- x[3]
# Create a vector for genotypic frequencies
Genotype_Freq <-vector("numeric",6)  

# Use Hardy-Weinberg equation to calculate genotypic frequencies from allelic frequencies
Genotype_Freq[1] <- p^2
Genotype_Freq[2] <- 2*p*q
Genotype_Freq[3] <- 2*p*r
Genotype_Freq[4] <- 2*q*r
Genotype_Freq[5] <- q^2
Genotype_Freq[6] <- r^2

# Print allelic frequencies
cat("Observed allelic frequencies:", "\n", "f(J) = ", p, "\n", "f(K) = ", q, "\n","f(L) = ", r, "\n")
cat("\n")

#i Print expected Hardy-Weinberg genotypic frequencies
cat("Expected Hardy-Weinberg genotypic frequencies:", "\n", "H-W f(JJ) = ", Genotype_Freq[1], "\n", "H-W f(JK) = ", Genotype_Freq[2], "\n", "H-W f(JL) = ", Genotype_Freq[3], "\n",  "H-W f(KL) = ", Genotype_Freq[4], "\n", "H-W f(KK) = ", Genotype_Freq[5], "\n", "H-W f(LL) = ", Genotype_Freq[6], "\n")
  cat("\n")

  # Return the output vector
  return(Genotype_Freq)
}

## ----tidy=TRUE-----------------------------------------------------------
AlleleFreq_2A(x=c(AA=75, AB=20, BB=100))

## ----tidy=TRUE-----------------------------------------------------------
HardyWeinberg_2A(x=c(p=0.4358974, q=0.5641026))

## ----tidy=TRUE-----------------------------------------------------------
AlleleFreq_3A(x=c(JJ=10,JK=11,JL=0,KL=9,KK=2,LL=22))

## ----tidy=TRUE-----------------------------------------------------------
HardyWeinberg_3A(x=c(p=0.287037, q=0.2222222, r=0.4907407))

## ----tidy=TRUE-----------------------------------------------------------
# Chaining functions together to get allelic frquencies and Hardy-Weinberg expected genotypic frequencies:

HardyWeinberg_2A(AlleleFreq_2A(x=c(AA=75, AB=20, BB=100)))

HardyWeinberg_3A(AlleleFreq_3A(x=c(JJ=10,JK=11,JL=0,KL=9,KK=2,LL=22)))

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate the increase in the frequency of a mutant allele through time
# NJG
# 21 November 2015

Mutation <- function(qo=0.5,u=0.000001,t=1:10) { 
  qt =  1 - (1 - qo)*exp(-u*t)
return(qt)
}

## ----tidy=TRUE-----------------------------------------------------------
Mutation(qo=0.5)


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate the change in allele frequency from migration
# NJG
# 21 November 2015

Migration <- function(p0=0.5, pm=0.9, m=0.1, t=1:10){
  pt <- (1 - m)^t * (p0 - pm) + pm
  return(pt)
}


## ----tidy=TRUE-----------------------------------------------------------
Migration(p0=0.1)

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate the change in allele frequency from inbreeding
# NJG
# 21 November 2015

Inbreeding <- function(p=0.3, F = 0.5){
  genotypes <- vector("numeric",3)
  q <- 1 - p
  genotypes[1] <- p^2*(1 - F) + p*F
  genotypes[2] <- 2*p*q*(1- F)
  genotypes[3] <- q^2*(1 - F) + q*F
  
  return(genotypes)
}


## ----tidy=TRUE-----------------------------------------------------------
Inbreeding()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate effective population size with a bottleneck
# NJG
# 21 November 2015

Bottleneck <- function(N=1:5){
  Ne <- 1/((1/length(N))*(sum(1/N)))
  return(Ne)
}


## ----tidy=TRUE-----------------------------------------------------------
Bottleneck()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate effective population size with a skewed sex ratio
# NJG
# 21 November 2015

SexRatio <- function(m=10, f=12){ 
  Ne <- (4*m*f)/(m + f)
  return(Ne)
}

## ----tidy=TRUE-----------------------------------------------------------
SexRatio()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate effective population size with limited dispersal
# NJG
# 21 November 2015

NatalDispersal <- function(d=10, x=1){
  Ne <- 4*pi*d*x
  return(Ne)
}


## ----tidy=TRUE-----------------------------------------------------------
NatalDispersal()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate probability of at least one occurrence with individual probability p and number of trials n
# NJG
# 21 November 2015

CompoundProb <- function(p=0.01, n=52){
  Prob <- 1 - (1 - p)^n
  return(Prob)
}

## ------------------------------------------------------------------------
CompoundProb()

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to Calculate and Print 7 steps of Natural Selection Equations
# NJG
# 21 November 2015

SevenSteps <- function(gen=c(50,50,100), w=c(0.4, 0.35, 0.2)){
  
# Step 1: Given Initial Genotype Counts AND Relative Fitness
  w <- w/max(w)
cat("Step 1: Given Initial Genotype Counts AND Relative Fitness","\n")
  cat(" #(AA) = ",gen[1], " #(AB) = ",gen[2]," #(BB) = ",gen[3],"\n")
  cat(" w1 = ",w[1], " w2 = ",w[2]," w3 = ",w[3],"\n")
  cat("\n")
  
# Step 2: Calculate Initial Genotype And Allelic Frequencies (p0, q0)
  gen <- gen/sum(gen)
cat("Step 2: Calculate Initial Genotype And Allelic Frequencies (p0, q0)", "\n")
  cat(" f(AA) = ",gen[1], " f(AB) = ",gen[2]," f(BB) = ",gen[3],"\n")
  p0 = gen[1] + 0.5*gen[2]
  q0 = gen[3] + 0.5*gen[2]
  cat(" f(A) =  ",p0, " f(B)  = ",q0,"\n")
  cat("\n")
  
# Step 3: Calculate Genotype Frequencies AFTER Random Mating
gen[1] <- p0^2
gen[2] <- 2*p0*q0
gen[3] <- q0^2

cat("Step 3: Calculate Genotype Frequencies AFTER Random Mating", "\n")
cat(" f(AA) = ",gen[1], " f(AB) = ",gen[2]," f(BB) = ",gen[3],"\n")
cat("\n")
  
# Step 4: Calculate Genotype Frequencies AFTER Selection

gen <- gen*w

wbar <- sum(gen)

cat("Step 4: Calculate Genotype Frequencies AFTER Selection", "\n")
cat(" f(AA) = ",gen[1], " f(AB) = ",gen[2]," f(BB) = ",gen[3],"\n")
cat("Mean fitness = ", wbar, "\n")
cat("\n")
  
# Step 5: Normalize Genotype Frequencies
gen <- gen/wbar
cat("Step 5: Normalize Genotype Frequencies", "\n")
cat(" f(AA) = ",gen[1], " f(AB) = ",gen[2]," f(BB) = ",gen[3],"\n")
cat("\n")

# Step 6: Calculate New Allelic Frequencies
  p1 = gen[1] + 0.5*gen[2]
  q1 = gen[3] + 0.5*gen[2]
cat("Step 6: Calculate New Allelic Frequencies", "\n")
  cat(" f(A) =  ",p1, " f(B)  = ",q1,"\n")
  cat("\n")
  
# Step 7: Calculate New Genotype Frequencies AFTER Random Mating
gen[1] <- p1^2
gen[2] <- 2*p1*q1
gen[3] <- q1^2

cat("Step 7: Calculate New Genotype Frequencies AFTER Random Mating", "\n")
cat(" f(AA) = ",gen[1], " f(AB) = ",gen[2]," f(BB) = ",gen[3],"\n")
cat("\n")

}

## ----tidy=TRUE-----------------------------------------------------------
SevenSteps()

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION Fisher engine to calculate changes in allelic frequency with selection in each generation
# NJG
# 21 November 2015
FisherEngine <- function(t=20,p0=0.1,w=c(1,1,0.5)){
  
# Create vectors for storing pi, wbar, and the 3 genotypes
  pvec <- vector(mode="numeric", length=(t + 1))
  wbar <- vector(mode="numeric", length=(t))
  gen <- vector(mode="numeric", length=3)
  
# Loop through the selection random mating calculations
  pvec[1] <- p0
  for (i in 2:(t + 1)){
    gen[1] <- pvec[i-1]^2*w[1]
    gen[2] <- 2*(1 - pvec[i - 1])*pvec[i - 1]*w[2]
    gen[3] <- (1 - pvec[i - 1])^2*w[3]
    
    wbar[i-1] <- sum(gen)
    gen <- gen/wbar[i-1]
    pvec[i] <- gen[1] + 0.5*gen[2]
  }
  # Graph p and wbar as a function of time
  par(mfrow=c(1,2))
  plot(x=1:(t+1),y=pvec,xlab="Generation",ylab="p",type="l", ylim=c(0,1),las=1)
  grid()
  plot(x=1:(t),y=wbar,xlab="Generation",ylab="Mean Fitness",type="l", ylim=c(0,1), las=1)
  grid()
  return(list(pvec,wbar))
  
}

## ----tidy=TRUE-----------------------------------------------------------
FisherEngine()

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION For Calculating Heritability From A Selective Breeding Experiment
# 22 November
# NJG

Heritability <- function(x=10,y=20,z=11){
  SelectionDifferential <- y - x
  cat("Selection Differential = ",SelectionDifferential,"\n")
  cat("\n")
  
  ResponseToSelection <- z - x
  cat("Response To Selection = ",ResponseToSelection, "\n")
  cat("\n")
  
  h2 <- ResponseToSelection/SelectionDifferential
  cat("Heritability = ",h2, "\n")
  cat("\n")
}

## ----tidy=TRUE-----------------------------------------------------------
Heritability()

## ----tidy=TRUE-----------------------------------------------------------

# FUNCTION for plotting results of common garden experiments
# 21 November 2015
# NJG

CommonGarden <- function(TraitMeans=c(40,20,70)){
  Pop1 <- rnorm(15,TraitMeans[1],10)
  Pop2 <- rnorm(15,TraitMeans[2],10)
  Pop3 <- rnorm(15,TraitMeans[3],10)
PopData <- c(Pop1,Pop2,Pop3)
Treatment <- rep(c("P1","P2","P3"),each=15)
par(mar=c(6,6,4,2))
BoxPlotDataFrame <- data.frame(x=Treatment,y=PopData)
boxplot(y~x,data=BoxPlotDataFrame,ylab="Trait Value",sub="Donor Source",cex.sub=2,cex.axis=1.5, cex.lab=2,col="bisque")
  }

## ----tidy=TRUE-----------------------------------------------------------
CommonGarden()

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION RecipPlot generates a plot of reciprocal transplant data
# NJG
# 21 November 2015
# The user specifies the mean trait values for each of the 4 treatments
# NJG
# 21 November 2015

RecipPlot <- function(ydata=c(10,20,30,5)){
  # set margins and set up axis locations
  par(mar=c(9,4,4,2))
  xdata <- c(1,1,2,2)
  
  # set up an empty plot
  plot(x=xdata,y=ydata,xlim=c(0.5,2.5),ylim=c(min(ydata)-5,max(ydata)+5),ann=F,axes=F,type="n")
  grid(nx=0,ny=10)
  # add x axis and labels
  axis(side=1,labels=c("Cold","Warm"),at=c(1,2),tick=T, cex.axis=1.5)
  mtext("Donor Population",side=1,cex=2,line=4)
  mtext("Trait Value",side=2,cex=2, line=1)
  box()
  
  # add lines and points
  lines(c(1,2),ydata[c(1,2)])
  lines(c(1,2),ydata[c(3,4)])
  points(c(1,2),ydata[c(1,2)],cex=4,pch=21,bg="skyblue")
  points(c(1,2),ydata[c(3,4)],cex=3,pch=22,bg="salmon")
  
  # add legend
  legend("topright",legend=c("Cold Recipient Site","Warm Recipient Site"),pch=c(21,22),pt.cex=2,pt.bg=c("skyblue","salmon"))
  
}
RecipPlot()


