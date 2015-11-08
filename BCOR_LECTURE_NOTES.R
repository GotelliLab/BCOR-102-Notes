## ----tidy=TRUE-----------------------------------------------------------
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

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate observed allele frequencies or a single gene with 2 alleles  

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

Mutation <- function(qo=0.5,u=0.000001,t=1:10) { 
  qt =  1 - (1 - qo)*exp(-u*t)
return(qt)
}

## ----tidy=TRUE-----------------------------------------------------------
Mutation(qo=0.5)


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate the change in allele frequency from migration

Migration <- function(p0=0.5, pm=0.9, m=0.1, t=1:10){
  pt <- (1 - m)^t * (p0 - pm) + pm
  return(pt)
}


## ----tidy=TRUE-----------------------------------------------------------
Migration(p0=0.1)

## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate the change in allele frequency from inbreeding

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

Bottleneck <- function(N=1:5){
  Ne <- 1/((1/length(N))*(sum(1/N)))
  return(Ne)
}


## ----tidy=TRUE-----------------------------------------------------------
Bottleneck()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate effective population size with a skewed sex ratio

SexRatio <- function(m=10, f=12){ 
  Ne <- (4*m*f)/(m + f)
  return(Ne)
}

## ----tidy=TRUE-----------------------------------------------------------
SexRatio()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate effective population size with limited dispersal

NatalDispersal <- function(d=10, x=1){
  Ne <- 4*pi*d*x
  return(Ne)
}


## ----tidy=TRUE-----------------------------------------------------------
NatalDispersal()


## ----tidy=TRUE-----------------------------------------------------------
# FUNCTION to calculate probability of at least one occurrence with individual probability p and number of trials n

CompoundProb <- function(p=0.01, n=52){
  Prob <- 1 - (1 - p)^n
  return(Prob)
}

## ------------------------------------------------------------------------
CompoundProb()

