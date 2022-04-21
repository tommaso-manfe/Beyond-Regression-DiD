library(tidyverse)
library(hrbrthemes) # theme ggplot

repository="C:/Users/tommy/OneDrive/Desktop/tesi/DID simulation/MC simulation"
setwd(repository)

##############################################################################
#DATASET
###############################################################################

# Sample size
n <- 10000
# pscore index (strength of common support)
Xsi.ps <- .75
# Proportion in each period
lambda <- 0.5
# NUmber of bootstrapped draws
nboot <- 199
#-----------------------------------------------------------------------------
# Mean and Std deviation of Z's without truncation
mean.z1 <- exp(0.25/2)
sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
mean.z2 <- 10
sd.z2 <- 0.54164
mean.z3 <- 0.21887
sd.z3 <-   0.04453
mean.z4 <- 402
sd.z4 <-  56.63891
#-----------------------------------------------------------------------------
# Gen covariates
x1 <- stats::rnorm(n, mean = 0, sd = 1)
x2 <- stats::rnorm(n, mean = 0, sd = 1)
x3 <- stats::rnorm(n, mean = 0, sd = 1)
x4 <- stats::rnorm(n, mean = 0, sd = 1)

z1 <- exp(x1/2)
z2 <- x2/(1 + exp(x1)) + 10
z3 <- (x1 * x3/25 + 0.6)^3
z4 <- (x1 + x4 + 20)^2

z1 <- (z1 - mean.z1)/sd.z1
z2 <- (z2 - mean.z2)/sd.z2
z3 <- (z3 - mean.z3)/sd.z3
z4 <- (z4 - mean.z4)/sd.z4

x <- cbind(x1, x2, x3, x4)
z <- cbind(z1, z2, z3, z4)
#-----------------------------------------------------------------------------
# Gen treatment groups
# Propensity score
pi <- 0.5
d  <- as.numeric(runif(n) <= pi)
#-----------------------------------------------------------------------------
# Generate aux indexes for the potential outcomes
index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)

# Create heterogenenous effects for the ATT, which is set approximately equal to zero
index.unobs.het <- d * (index.lin)
index.att <- -10*x1+10*x2-10*x3-10*x4
index.att=index.att-mean(index.att[d==1])
ATT[i] <- mean(index.att[d==1])

#This is the key for consistency of outcome regression
index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
#v is the unobserved heterogeneity
v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)

#Gen realized outcome at time 0
y0 <- index.lin + v + stats::rnorm(n)

# gen outcomes at time 1
# First let's generate potential outcomes: y_1_potential
y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
  index.trend #this is for the trend based on X

y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
  index.trend + #this is for the trend based on X
  index.att # This is the treatment effects

# Gen realized outcome at time 1
y1 <- d * y11 + (1 - d) * y10

# Generate "T"
# Generate "T"
post <- as.numeric(stats::runif(n) <= lambda)
#ti <- stats::plogis(- x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4)
#post  <- as.numeric(runif(n) <= ti)
# observed outcome
y <- post * y1 + (1 - post) * y0
#-----------------------------------------------------------------------------
#Gen id
id <- 1:n
#-----------------------------------------------------------------------------
# Put in a long data frame
dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                x1 = z1, x2= z2, x3 = z3, x4 = z4))
dta_long <- dta_long[order(dta_long$id),]

################################################################################
datax2=as.data.frame(x2[d==1 & post==0])
colnames(datax2)=c('x2')
datax2$post=0
temp=as.data.frame(x2[d==1 & post==1])
colnames(temp)=c('x2')
temp$post=1
datax2=rbind(datax2,temp)

# plot multilevel density plot
# plot multileve density with facet_warp
treated <- ggplot(data=datax2, aes(x=x2, group=post, fill=post)) +
  geom_density(adjust=1.5, alpha=.4)+
  theme_ipsum()+
  scale_x_continuous(limits = c(-3, 3))+
  geom_vline(data = datax2, aes(xintercept = mean(x2[post==1]), 
                                color = post), size=1.5)+
  geom_vline(data = datax2, aes(xintercept = mean(x2[post==0]) 
  ), size=1.5)+
ggsave("EXP0_treated.png")

datax2=as.data.frame(x2[d==0 & post==0])
colnames(datax2)=c('x2')
datax2$post=0
temp=as.data.frame(x2[d==0 & post==1])
colnames(temp)=c('x2')
temp$post=1
datax2=rbind(datax2,temp)

# plot multilevel density plot
# plot multileve density with facet_warp
controls <- ggplot(data=datax2, aes(x=x2, group=post, fill=post)) +
  geom_density(adjust=1.5, alpha=.4)+
  theme_ipsum()+
  scale_x_continuous(limits = c(-3, 3))+
  geom_vline(data = datax2, aes(xintercept = mean(x2[post==1]), 
                                color = post), size=1.5)+
  geom_vline(data = datax2, aes(xintercept = mean(x2[post==0]) 
  ), size=1.5)

ggsave("EXP0_controls.png")








