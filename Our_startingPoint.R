data("faithful")
head(faithful)
# Point 1 Creating the dataset ---------------------------------------------------------
handmade.em <- function(y,plot_flag = T){ #plot_flag = T -> make a plot
  k = 5
  waiting_max = max(faithful$waiting)
  waiting_min = min(faithful$waiting)
  mu = floor(runif(k,waiting_min,waiting_max))
  sigma = floor(runif(k,1,sqrt(var(faithful$waiting)))) #to evaluate
  p = c(0.2,0.2,0.2,0.2,0.2) #1/k
  
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
    # E step - 'get the responsabilities' - optimal E-step distribution
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
  
    # M step - '' 
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
  
  res <- data.frame(res)
  names(res) <- c("iteration","p1","p2","p3","p4","p5","mu1","mu2","mu3","mu4","mu5","sigma1","sigma2","sigma3","sigma4", "sigma5","deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}

dataset <- handmade.em(faithful$waiting) #taking only the 'waiting' column
                       #initial values for our parameters - ALL RANDOM
#show our dataset
round(dataset$parameters, 3)
hem_fit$deviance



# Point 2 get the sampling ------------------------------------------------
#since we can apply the LLN to any sample with size bigger than 40, we chose the choose reasonably small and the second one bigger
#(to evaluate that values)
n1 = 15
n2 = 100


