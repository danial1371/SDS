data("faithful")
?faithful
head(faithful)
# Handmade EM4MoG ---------------------------------------------------------
#1 dim case

handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T){ #plot_flag = T -> make a plot
  # Init / 2 components only
  cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
  
  # ------ THESE ARE THE STARTING POINT OF OUR EVALUATION ----------------
  #the likelihood is build as sum of 2 components (it's a vector)
  like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2]) #evaluation of our likelihood (for 2 components)
  deviance <- -2*sum(log(like)) #applying the definition of the deviance
  res      <- matrix(NA,n_iter + 1, 8) #initialization of the matrix result (n_iter+1,8)
  #The matrix has 8 cols because we have 4 paramters with 2 components 
 
  res[1,]  <- c(0, p, mu, sigma, deviance)  #0= starting iterations; p,mu,sigma=initial values, deviance = value of the deviance
  
  for (iter in 1:n_iter) {
    # E step - 'get the responsabilities' - optimal E-step distribution
    # FOR EACH OBSERVATION
    d1 <- p[1]*dnorm(y, mu[1], sigma[1])
    d2 <- p[2]*dnorm(y, mu[2], sigma[2])
    #normalize it
    r1 <- d1/(d1 + d2) 
    r2 <- 1 - r1
    
    # M step - '' 
    # FOR EACH OBSERVATION
    p[1]     <- mean(r1)
    mu[1]    <- sum(r1*y)/sum(r1) #weighted mean of the observations
    sigma[1] <-sqrt( sum(r1*(y^2))/sum(r1) - (mu[1])^2 ) #weighted variance for the std
    
    p[2]     <- 1 - p[1]
    mu[2]    <- sum((r2)*y)/sum((r2))
    sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    
    #-----------UPDATE THE TARGETS------------------
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot the distributions
    if (plot_flag){
      hist(y, prob = T, breaks = 30, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             #classes
             col = cols[ (dnorm(y,mu[1],sigma[1]) > dnorm(y,mu[2],sigma[2])) + 1])
      curve(p[1]*dnorm(x,mu[1],sigma[1]) + p[2]*dnorm(x,mu[2],sigma[2]),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  names(res) <- c("iteration","p1","p2","mu1","mu2","sigma1","sigma2","deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}

                      #y = data;
hem_fit <- handmade.em(faithful$waiting, #taking only the 'waiting' column
                       #initial values for our parameters - ALL RANDOM
#if we want to choose them, look at summary(faithful)
                       p      = c(.5,.5), #p= proportions
                       mu     = c(45,55), #mu= centers -> ARBITRARY
                       sigma  = c(8,8), #sigma= std
                       n_iter = 20) #n_iter = number of iterations
round(hem_fit$parameters, 3 )
hem_fit$deviance


# Handmade EM: Different Initial Values -----------------------------------

require(ggplot2)
require(gridExtra)
plotConvMC <- function(df, title = NULL){
  G      <- (ncol(df)-2)/3
  df$rep <- as.factor(df$rep)
  graf   <- vector("list", ncol(df)-2)
  for (j in (2:(ncol(df)-1))) {
    grafj <- ggplot(df) + geom_line(aes_string(df[,1],df[,j], 
                                               color = df[,ncol(df)])) +
      xlab("iteration") + ylab(names(df[j])) + theme(legend.position = "none")
    graf[[j-1]] <- grafj
  }
  do.call("grid.arrange", c(graf, ncol = 3, top = title))
}

D.em <- NULL
set.seed(1234)
for (m in (1:10)) {
  p1 <- runif(1,0.1,0.9)
  df.em <- handmade.em(faithful$waiting, 
                       p  = c(p1, 1 - p1), 
                       mu = rnorm(2, 70, 15), sigma = rlnorm(2, 2, 0.7),
                       n_iter = 50, plot_flag = F)$res
  df.em$rep <- m
  D.em <- rbind(D.em,df.em)
}
plotConvMC(D.em)

# Up to some permutation (labels switching), all the runs converge to the same solution.
# A very poor initial guess may lead to a very poor convergence of EM.

handmade.em(faithful$waiting, p = c(0.5, 0.5), mu = c(30, 40), sigma = c(2, 10), n_iter = 5)

# In this example, one of the variances converges essentially to 0... so be carfull while choosing the random initialization!

