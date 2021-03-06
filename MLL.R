data("faithful")
head(faithful)
# Point 1 Creating the dataset ---------------------------------------------------------
handmade.em <- function(y,plot_flag = T){ #plot_flag = T -> make a plot
  k = 5
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
  n_iter = 50  #n_iter = number of iterations
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
    mu[3]    <- sum(r3*y)/sum(r3) #weighted mean of the observations
    sigma[3] <-sqrt(sum(r3*(y^2))/sum(r3) - (mu[3])^2 ) #weighted variance for the std
    
    p[4]     <- mean(r4)
    mu[4]    <- sum((r4)*y)/sum((r4))
    sigma[4] <- sqrt(sum(r4*(y^2))/sum(r4) - (mu[4])^2)

    p[5]     <- mean(r5)
    mu[5]    <- sum((r5)*y)/sum((r5))
    sigma[5] <- sqrt(sum(r5*(y^2))/sum(r5) - (mu[5])^2)
    
    #normalize the vector p to get the summation equal to 1
    p = p/sum(p)
    
    #-----------UPDATE THE TARGETS------------------
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2]) + p[3]*dnorm(y, mu[3], sigma[3])+ p[4]*dnorm(y, mu[4], sigma[4])+ p[5]*dnorm(y, mu[5], sigma[5])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
  }
  
 # res <- data.frame(res)
 # names(res) <- c("iteration","p1","p2","p3","p4","p5","mu1","mu2","mu3","mu4","mu5","sigma1","sigma2","sigma3","sigma4", "sigma5","deviance")
 # out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  #return(out)
  return(res)
}

dataset <- handmade.em(faithful$waiting) #taking only the 'waiting' column
                       #initial values for our parameters - ALL RANDOM
#show our dataset
#round(dataset$parameters, 3)
#dataset$deviance



# Point 2 get the sampling ------------------------------------------------
#since we can apply the LLN to any sample with size bigger than 40, we chose the choose reasonably small and the second one bigger
#(to evaluate that values)
n1 = 18
n2 = 100
guassian1 = rnorm(n1,0,1) # firat one is standard normal distribution 
guassian2 = rnorm(n1,dataset[50,7],dataset[50,12]) # second argument is last mu1 and third argument is last sigma1 and first is sample size (n1)
guassian3 = rnorm(n1,dataset[50,8],dataset[50,13])
guassian4 = rnorm(n1,dataset[50,9],dataset[50,14])
guassian5 = rnorm(n1,dataset[50,10],dataset[50,15])
guassian6 = rnorm(n1,dataset[50,11],dataset[50,16])

n1_prob_density = 0.5*(guassian1)+(0.1)*(guassian2)+(0.1)*(guassian3)+(0.1)*(guassian4)+(0.1)*(guassian5)+(0.1)*(guassian6)
guassian1_n2 = rnorm(n2,0,1) # firat one is standard normal distribution 
guassian2_n2 = rnorm(n2,dataset[50,7],dataset[50,12]) # second argument is last mu1 and third argument is last sigma1 and first is sample size (n1)
guassian3_n2 = rnorm(n2,dataset[50,8],dataset[50,13])
guassian4_n2 = rnorm(n2,dataset[50,9],dataset[50,14])
guassian5_n2 = rnorm(n2,dataset[50,10],dataset[50,15])
guassian6_n2= rnorm(n2,dataset[50,11],dataset[50,16])

n2_prob_density = 0.5*(guassian1_n2)+(0.1)*(guassian2_n2)+(0.1)*(guassian3_n2)+(0.1)*(guassian4_n2)+(0.1)*(guassian5_n2)+(0.1)*(guassian6_n2)

# Point 3 Simulating and AIC ------------------------------------------------
n1 = 18 

suppressMessages(require(mixtools, quietly = T)) # Package and function ?rnormmix 

# the first parameter is lambda and we use our p normalized it to 0.5 in order to obtain at summation 1 and keep 0.5 standard weight
# for normal standard at out main formula of bart density 
# for mean and variance we use the first normal as standard normal and for the rest we replace the mean of all results drive from MLE 
XX1 <- rnormmix(n1,lambda = c(0.5, dataset[50,2:6]/2), mu = c(0, colMeans(dataset[,7:11])) , sigma = c(1, colMeans(dataset[,12:16]) )) 
# first bart just include normal studard and k=5 we make summation on normal standard with weigth equal to 0.5 and for the rest 0.5/4  
XX2 <- rnormmix(n1,lambda = c(0.5, dataset[50,3:6]/2), mu = c(0,colMeans(dataset[,8:11])) , sigma = c(1, colMeans(dataset[,13:16]) ))
XX3 <- rnormmix(n1,lambda = c(0.5, dataset[50,4:6]/2), mu = c(0, colMeans(dataset[,9:11])) , sigma = c(1, colMeans(dataset[,14:16]) ))
XX4 <- rnormmix(n1,lambda = c(0.5, dataset[50,5:6]/2), mu = c(0 ,colMeans(dataset[,10:11])), sigma = c(1, colMeans(dataset[,15:16]) ))
XX5 <- rnormmix(n1,lambda = c(0.5, dataset[50,6]/2)  , mu = c(0 ,mean(dataset[,11]))   , sigma = c(1, mean(dataset[,16]) ))


# Point 4 Applying MLE to each simulated data  ------------------------------------------------
library(KScorrect)
library(bbmle)
# LogLikelihood -----------------------------------------------------------
LL = function(par0,data){
  #R = dmixnorm(x = data, mean = par0[1:6] , sd = par0[7:12],pro = par0[13:18])
  x = data
  R = par0[13]*dnorm(x, 0, 1) +
    par0[14]*dnorm(x,par0[2], par0[8]) + 
    par0[15]*dnorm(x, par0[3], par0[9]) +
    par0[16]*dnorm(x, par0[4], par0[10]) + 
    par0[17]*dnorm(x, par0[5],par0[11]) +
    par0[18]*dnorm(x, par0[6], par0[12])
  -sum(log(R)) #we need the min value
}
n2 = 10
mean = c(0, colMeans(dataset[,7:11]))
pro = c(0.5, dataset[50,2:6]/2)
sd = c(1, colMeans(dataset[,12:16]) )
XX1 = rmixnorm(n2, mean = c(0, colMeans(dataset[,7:11])),pro = c(0.5, dataset[50,2:6]/2) , sd = c(1, colMeans(dataset[,12:16]) ))

par0 = c(mean =  c(0, colMeans(dataset[,7:11])), sd = c(1, colMeans(dataset[,12:16]) ),pro = c(0.5, dataset[50,2:6]/2))
fit = optim(par = par0, fn = LL, data = XX1,method ="L-BFGS-B",lower = c(rep(0,6),rep(0,12)) , upper = c(rep(10000,12),rep(1,6)) )

# SIMUALTIONS -------------------------------------------------------------

M=50 #numb of simulations
XX <- rep(NA, M)
simulations <- rep(NA, M)
max = 18*M
for (i in 1:M){
  XX = rmixnorm(n2, mean = c(0, colMeans(dataset[,7:11])),pro = c(0.5, dataset[50,2:6]/2) , sd = c(1, colMeans(dataset[,12:16]) ))
  for (j in 1:max){
    fit = optim(par = par0, fn = LL, data = XX,method ="L-BFGS-B",lower = c(rep(0,6),rep(0,12)) , upper = c(rep(10000,12),rep(1,6)))
    simulations[j] = fit$par
  }
}






