data("faithful")
head(faithful)
options(digits = 10) 
# Point 1 Creating the dataset ---------------------------------------------------------
handmade.em <- function(y,plot_flag = T){ #plot_flag = T -> make a plot
  k = 5
  n_iter = 10000  # number of iterations
  waiting_max = max(faithful$waiting)
  waiting_min = min(faithful$waiting)
  lower_quantile =  quantile(faithful$waiting, 0.3) 
  upper_quantile =  quantile(faithful$waiting, 0.7)
  mu = floor(runif(k,lower_quantile,upper_quantile))
  square_root = sqrt(var(faithful$waiting))
  sigma = floor(runif(k,square_root-3,square_root +3  )) #to evaluate
  p = c(rep(1/k,5)) #1/k
  
  # Init / 2 components only
  cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
  # ------ THESE ARE THE STARTING POINT OF OUR EVALUATION ----------------
  #the likelihood is build as sum of 2 components (it's a vector)
  like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2]) + p[3]*dnorm(y, mu[3], sigma[3])+ p[4]*dnorm(y, mu[4], sigma[4])+ p[5]*dnorm(y, mu[5], sigma[5])#evaluation of our likelihood (for 2 components)
  
  deviance <- -2*sum(log(like)) #applying the definition of the deviance
  res      <- matrix(NA,n_iter + 1, 17) #initialization of the matrix result (n_iter+1,20)
  #The matrix has 20 cols because we have 4 parameters with 5 components 
  
  res[1,]  <- c(0, p, mu, sigma, deviance)  #0= starting iterations; p,mu,sigma=initial values, deviance = value of the deviance
  
  for (iter in 1:n_iter) {
    # E step - Compute the expectations, or responsibilities, of the latent variables
    # FOR EACH OBSERVATION
    d1 <- p[1]*dnorm(y, mu[1], sigma[1])
    d2 <- p[2]*dnorm(y, mu[2], sigma[2])
    d3 <- p[3]*dnorm(y, mu[3], sigma[3])
    d4 <- p[4]*dnorm(y, mu[4], sigma[4])
    d5 <- p[5]*dnorm(y, mu[5], sigma[5])
    
    #normalize it
    r1 <- d1/(d1 + d2 + d3 + d4 + d5) 
    r2 <- d2/(d1 + d2 + d3 + d4 + d5) 
    r3 <- d3/(d1 + d2 + d3 + d4 + d5) 
    r4 <- d4/(d1 + d2 + d3 + d4 + d5) 
    r5 <- d5/(d1 + d2 + d3 + d4 + d5) 
  
    # M step - Compute the maximum likelihood parameters given these responsibilities 
    # FOR EACH OBSERVATION
    p[1]     <- mean(r1)
    mu[1]    <- sum(r1*y)/sum(r1) #weighted mean of the observations
    sigma[1] <-sqrt( sum(r1*(y^2))/sum(r1) - (mu[1])^2 ) #weighted variance for the std
    
    p[2]     <- mean(r2)
    mu[2]    <- sum((r2)*y)/sum((r2))
    sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    
    p[3]     <- mean(r3)
    mu[3]    <- sum(r3*y)/sum(r3) 
    sigma[3] <-sqrt(sum(r3*(y^2))/sum(r3) - (mu[3])^2 ) 
    
    p[4]     <- mean(r4)
    mu[4]    <- sum((r4)*y)/sum((r4))
    sigma[4] <- sqrt(sum(r4*(y^2))/sum(r4) - (mu[4])^2)

    p[5]     <- mean(r5)
    mu[5]    <- sum((r5)*y)/sum((r5))
    sigma[5] <- sqrt(sum(r5*(y^2))/sum(r5) - (mu[5])^2)
    
    #normalize the vector p to get the summation equal to 1
    p = p/sum(p)
    
    #-----------UPDATE THE TARGETS------------------
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2]) + p[3]*dnorm(y, mu[3], sigma[3])+ p[4]*dnorm(y, mu[4], sigma[4])+ p[5]*dnorm(y, mu[5], sigma[5])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
   
    #check the convergence of the loglikelihood
    # dev_diff = abs((res[iter,17] - res[iter-1,17]))
    # print(dev_diff)
    # if(dev_diff < 1e-5) {
    #   #return(res)
    #   break
    # }
  }
  return(res)
}

result <- handmade.em(faithful$waiting) #taking only the 'waiting' column
                       #initial values for our parameters - ALL RANDOM


#Change the matrix result into a dataframe
dataset <- data.frame(result)
names(dataset) <- c("iteration","p1","p2","p3","p4","p5","mu1","mu2","mu3","mu4","mu5","sigma1","sigma2","sigma3","sigma4", "sigma5","deviance")

dataset$deviance

# Point 2 get the sampling ------------------------------------------------
#since we can apply the LLN to any sample with size bigger than 40, we chose the first reasonably small and the second one bigger
#saving the last elements from our dataset
mu1= tail(dataset$mu1,1)
mu2= tail(dataset$mu2,1)
mu3= tail(dataset$mu3,1)
mu4= tail(dataset$mu4,1)
mu5= tail(dataset$mu5,1)

sigma1 = tail(dataset$sigma1,1)
sigma2 = tail(dataset$sigma2,1)
sigma3 = tail(dataset$sigma3,1)
sigma4 = tail(dataset$sigma4,1)
sigma5 = tail(dataset$sigma5,1)

p = c(0.5, tail(dataset$p1,1)/2,tail(dataset$p2,1)/2,tail(dataset$p3,1)/2,tail(dataset$p4,1)/2,tail(dataset$p5,1)/2)

#FIRST SAMPLE
n1 = 18
guassian1_n1 = rnorm(n1,0,1) # first one is standard normal distribution 
# second argument is last mu1 and third argument is last sigma1 and first is sample size (n1)
guassian2_n1 = rnorm(n1,mu1,sigma1) 
guassian3_n1 = rnorm(n1,mu2,sigma2)
guassian4_n1 = rnorm(n1,mu3,sigma3)
guassian5_n1 = rnorm(n1,mu4,sigma4)
guassian6_n1 = rnorm(n1,mu5,sigma5)

n1_prob_density = 0.5*(guassian1_n1)+
                  p[2]*(guassian2_n1)+
                  p[3]*(guassian3_n1)+
                  p[4]*(guassian4_n1)+
                  p[5]*(guassian5_n1)+
                  p[6]*(guassian6_n1)

# SECOND SAMPLE
n2 = 200
guassian1_n2 = rnorm(n2,0,1) # first one is standard normal distribution 
# second argument is last mu1 and third argument is last sigma1 and first is sample size (n1)
guassian2_n2 = rnorm(n2,mu1,sigma1) 
guassian3_n2 = rnorm(n2,mu2,sigma2)
guassian4_n2 = rnorm(n2,mu3,sigma3)
guassian5_n2 = rnorm(n2,mu4,sigma4)
guassian6_n2 = rnorm(n2,mu5,sigma5)

n2_prob_density = 0.5*(guassian1_n2)+
                  p[2]*(guassian2_n2)+
                  p[3]*(guassian3_n2)+
                  p[4]*(guassian4_n2)+
                  p[5]*(guassian5_n2)+
                  p[6]*(guassian6_n2)

#Simulating ------------------------------------------------
suppressMessages(require(mixtools, quietly = T)) # Package and function
# the first parameter is lambda and we use our p normalized it to 0.5 in order to obtain at summation 1 and keep 0.5 standard weight
# for mean and variance we use the first normal as standard normal and for the rest we replace the mean of all results drive from MLE 
XX1 <- rnormmix(n1,lambda = c(0.5, dataset[50,2:6]/2), mu = c(0, colMeans(dataset[,7:11])) , sigma = c(1, colMeans(dataset[,12:16]) )) 
# first bart just include normal standard and k=5 we make summation on normal standard with weight equal to 0.5 and for the rest 0.5/4  
XX2 <- rnormmix(n1,lambda = c(0.5, dataset[50,3:6]/2), mu = c(0,colMeans(dataset[,8:11])) , sigma = c(1, colMeans(dataset[,13:16]) ))
XX3 <- rnormmix(n1,lambda = c(0.5, dataset[50,4:6]/2), mu = c(0, colMeans(dataset[,9:11])) , sigma = c(1, colMeans(dataset[,14:16]) ))
XX4 <- rnormmix(n1,lambda = c(0.5, dataset[50,5:6]/2), mu = c(0 ,colMeans(dataset[,10:11])), sigma = c(1, colMeans(dataset[,15:16]) ))
XX5 <- rnormmix(n1,lambda = c(0.5, dataset[50,6]/2)  , mu = c(0 ,mean(dataset[,11]))   , sigma = c(1, mean(dataset[,16]) ))


#Applying MLE to the first simulated data (n1) ------------------------------------------------


# LogLikelihood -----------------------------------------------------------
#par0 is the vector of parameter values that the optimization algorithm wants to test
#data = the data for which the NLL is calculated

LL = function(par0,data){
  #R = dmixnorm(x = data, mean = par0[1:6] , sd = par0[7:12],pro = par0[13:18])
  x = data
  #the dnorm returns the probability density of the data assuming a Normal distribution with given mean and standard deviation
  # the log=TRUE  the logarithm of the probability density
  Likelihood =  par0[13]*dnorm(x, 0, 1) +
                par0[14]*dnorm(x,par0[2], par0[8]) + 
                par0[15]*dnorm(x, par0[3], par0[9]) +
                par0[16]*dnorm(x, par0[4], par0[10]) + 
                par0[17]*dnorm(x, par0[5],par0[11]) +
                par0[18]*dnorm(x, par0[6], par0[12])
  NLL = -sum(log(Likelihood)) #we need the min value
  return (NLL)
}

# One trial (DELETE WHEN THE M SIMULATIONS WORK) -------------------------------------------------------------
pro = c(0.5, dataset[50,2:6]/2)
mean = c(0, colMeans(dataset[,7:11]))
sd = c(1, colMeans(dataset[,12:16]) )
XX1 = rmixnorm(n1, mean = c(0, colMeans(dataset[,7:11])),pro = c(0.5, dataset[50,2:6]/2) , sd = c(1, colMeans(dataset[,12:16]) ))

par0 = c(mean =  c(0, colMeans(dataset[,7:11])), sd = c(1, colMeans(dataset[,12:16]) ),pro = c(0.5, dataset[50,2:6]/2))
#minimize the LL using optim() that needs the initial values for each parameter,the function calculating LL
#and the arguments that will be passed to the objective function
fit = optim(par = par0, fn = LL, data = XX1,method ="L-BFGS-B",lower = c(rep(0,6),rep(0,12)) , upper = c(rep(10000,12),rep(1,6)) )

#this is the only output that we want
fit$par #This shows the MLEs of the parameters.


# M SIMUALTIONS -------------------------------------------------------------

M=50 #number of simulations
#XX <- rep(NA, M)
simulations <- rep(NA, M)
max = 18*M
for (i in 1:M){
  XX = rmixnorm(n2, mean = c(0, colMeans(dataset[,7:11])),pro = c(0.5, dataset[50,2:6]/2) , sd = c(1, colMeans(dataset[,12:16]) ))
  for (j in 1:max){
    fit = optim(par = par0, fn = LL, data = XX,method ="L-BFGS-B",lower = c(rep(0,6),rep(0,12)) , upper = c(rep(10000,12),rep(1,6)))
    simulations[j] = fit$par
  }
}





# Point 3 get the simulations ---------------------------------------------
library(KScorrect)
library(bbmle)
#first simulations for the sample sizes n1 and n2
XX1 = rmixnorm(n1, mean = c(0,mu1,mu2,mu3,mu4,mu5),pro = p, sd = c(1,sigma1,sigma2,sigma3,sigma4,sigma5))
XX1_df = data.frame(XX1)
XX2 = rmixnorm(n2, mean = c(0,mu1,mu2,mu3,mu4,mu5),pro = p, sd = c(1,sigma1,sigma2,sigma3,sigma4,sigma5))
XX2_df = adata.frame(XX2)

#pick random values of size n1 and n2 from our imput dataset in order to fit the simulation
sample1 = as.matrix(faithful$waiting)[sample(nrow(as.matrix(faithful$waiting)), n1, replace = FALSE)]
sample1_df = data.frame(sample1)
sample2 = as.matrix(faithful$waiting)[sample(nrow(as.matrix(faithful$waiting)), n2, replace = FALSE)]
sample2_df = data.frame(sample2)

#WE DON'T KNOW IF WE HAVE TO USE THIS OR XX1
gauss1_n1 = data.frame(rnorm(n1,0,1)*0.5)
gauss2_n1 = data.frame(rnorm(n1,mu2,sigma2)*p[2])
gauss3_n1 = data.frame(rnorm(n1,mu3,sigma3)*p[3])
gauss4_n1 = data.frame(rnorm(n1,mu3,sigma3)*p[4])
gauss5_n1 = data.frame(rnorm(n1,mu4,sigma4)*p[5])
gauss6_n1 = data.frame(rnorm(n1,mu5,sigma5)*p[6])

#a. AIC ---------------------------------------------------------------------
#for the smaller sample size n1
#reg1 = lm(formula = sample1_df$sample1 ~ XX1_df$XX1)

reg1 = lm(formula = sample1_df$sample1  ~ gauss1_n1$guassian1_n1 + gauss2_n1$guassian2_n1 + gauss3_n1$guassian3_n1 + gauss4_n1$guassian4_n1 +gauss5_n1$guassian5_n1 + gauss6_n1$guassian6_n1 )
AIC(reg1)

M = 50

#for the smaller sample size n2
reg2 = lm(formula = sample2_df$sample2 ~ XX2_df$XX2)
AIC(reg2)

#b. BIC ---------------------------------------------------------------------
#for the smaller sample size n1
BIC(reg1)


#for the smaller sample size n2
BIC(reg2)




#c. SPLITTING 50-50 ---------------------------------------------------------------------
n = length(XX1)
data.points = 1:n
data.points = sample(data.points) # Permute randomly
train = data.points[1:floor(n/2)] # First random half is training
test = data.points[-(1:floor(n/2))] # 2nd random half is testing

candidate.component.numbers <- 2:10
loglikes <- vector(length=1+length(candidate.component.numbers))
# k=1 needs special handling
mu<-mean(XX1[train]) # MLE of mean
sigma <- sd(XX1[train])*sqrt((n-1)/n) # MLE of standard deviation
loglikes[1] <- sum(dnorm(XX1[test],mu,sigma,log=TRUE))
for (k in candidate.component.numbers) {
  mixture <- normalmixEM(XX1[train],k=k,maxit=400,epsilon=1e-2)
  loglikes[k] <- loglike.normalmix(XX1[test],mixture=mixture)
}



#d. SPLITTING 70-30 ---------------------------------------------------------------------







#e. SPLITTING 30-70 ---------------------------------------------------------------------






#f. 5-fold CROSS VALIDATION ---------------------------------------------------------------------






#g. 10-fold CROSS VALIDATION ---------------------------------------------------------------------





#h. Wass-score ---------------------------------------------------------------------






