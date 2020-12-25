# Homemade FWD sampling ---------------------------------------------------
#you sample your data (e.g the discrete latent variable for a fixed point of components) and one selecting the component, you generate your data

M   <- 100
# M   <- 1000
out <- rep(NA, M)
for (m in 1:M){
  #we are in the modeling stage so some information are known
  #select the label with a proportion you pick (the 0.4,0.6 are known according to the sample)
  lab = sample(c("Female", "Male"), size = 1, prob = c(.4,.6))
  if (lab == "Female") 
    #known information
    out[m] <- rnorm(n = 1, mean = 162, sd = sqrt(5))
  else 
    out[m] <- rnorm(n = 1, mean = 176, sd = sqrt(8))
}
hist(out, main = "", xlab = "")

head(out)
#now we dont' know each observation from which subpopulation come from; this is
#what we wanna re build from a clustering prospective

# Bart's density ----------------------------------------------------------
#IT USES DENSITY ESTIMATION AS BANCKMARK
########################################################################
# X ~ Bart Simpson's density (a.k.a. as "The Claw")   
      #6 components
#     f(x) = 0.5*dnorm(0,1) + sum(mu in 0:4){0.1*dnorm((mu/2)-1, 0.1)} #proportion
# 
# Mixture of 6 Normal densities with:
# pi = c(0.5,  0.1,  0.1,  0.1,  0.1,  0.1) #centers 
# mu = c(0.0, -1.0, -0.5,  0.0,  0.5,  1.0) #expectations 
# s  = c(1.0,  0.1,  0.1,  0.1,  0.1,  0.1) #standard dv 
########################################################################

# Some colors, just in case...
?rgb
colos <- c(rgb(32/255, 74/255, 135/255, 0.7),
           rgb(204/255, 0, 0, 0.7),
           rgb(200/255, 141/255, 0, 0.7),
           rgb(78/255, 153/255, 6/255, 0.7),
           rgb(32/255, 74/255, 135/255, 0.3),
           rgb(204/255, 0, 0, 0.3),
           rgb(200/255, 141/255, 0, 0.3),
           rgb(78/255, 153/255, 6/255, 0.3))

# Package and function
suppressMessages(require(mixtools, quietly = T))
?rnormmix

n <- 250
XX <- rnormmix(n, 
               lambda = c(0.5, rep(0.1,5)), 
               mu     = c(0, ((0:4)/2)-1), 
               sigma  = c(1, rep(0.1,5)) )

# Make an histogram of the data
hist(XX, prob = T, col = gray(.8), border = NA, xlab = "x",
     main = paste("Data from Bart's density",sep=""),
     sub = paste("n = ", n, sep = ""),
     breaks = 50)
# Show the data points
rug(XX, col = rgb(0,0,0,.5))

# Plot the true density
true.den = function(x) 0.5*dnorm(x, 0, 1) + 
  0.1*dnorm(x,-1.0, 0.1) + 0.1*dnorm(x, -0.5, 0.1) +
  0.1*dnorm(x, 0.0, 0.1) + 0.1*dnorm(x,  0.5, 0.1) +
  0.1*dnorm(x, 1.0, 0.1)
curve(true.den, col = rgb(1,0,0,0.4), lwd = 3, n = 500, add = TRUE)

# Kernel density estimate
lines(density(XX),            col = colos[3], lwd = 3)   # Oversmoothing
lines(density(XX, bw = .08),  col = colos[4], lwd = 3)   # Just Right
lines(density(XX, bw = .008), col = colos[5], lwd = 3)   # Undersmoothing

# Add a legend
legend("topright", c("True","Over", "Just right", "Under"), lwd = 5,
       col = c(rgb(1,0,0,0.4), colos[3], colos[4],colos[5]), cex = 0.8, bty = "n")

# Handmade EM4MoG ---------------------------------------------------------
#1 dim case
handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T){
  # Init / 2 components only
  cols     <- c(rgb(1,0,0,.3), rgb(0,1,0,.3))
  like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2]) #evaluation of our likelihood
  deviance <- -2*sum(log(like)) #applying the definition of the deviance
  res      <- matrix(NA,n_iter + 1, 8) #initialization of the matrix result (n_iter+1,8)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  for (iter in 1:n_iter) {
    # E step - 'get the responsabilities' - optimal E-step distribution
    d1 <- p[1]*dnorm(y, mu[1], sigma[1])
    d2 <- p[2]*dnorm(y, mu[2], sigma[2])
    r1 <- d1/(d1 + d2) #normalize it
    r2 <- 1 - r1
    
    # M step - ''
    p[1]     <- mean(r1)
    mu[1]    <- sum(r1*y)/sum(r1)
    sigma[1] <-sqrt( sum(r1*(y^2))/sum(r1) - (mu[1])^2 )
    p[2]     <- 1 - p[1]
    mu[2]    <- sum((r2)*y)/sum((r2))
    sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag){
      hist(y, prob = T, breaks = 30, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
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

data("faithful")
?faithful
hem_fit <- handmade.em(faithful$waiting, #taking only the 'waiting' column
                       #initial values for our parameters - ALL RANDOM
                       p      = c(.5,.5), 
                       mu     = c(45,55), 
                       sigma  = c(8,8), 
                       n_iter = 20)
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

# In this example, one of the variances converges essentially to 0...

# mixtools ----------------------------------------------------------------

?normalmixEM
#input data
mt_fit = normalmixEM(faithful$waiting)
summary(mt_fit)
names(mt_fit)

plot(mt_fit)
class(mt_fit)
?plot.mixEM
plot(mt_fit, density = TRUE, w = 1.1)

# Assign classes
(mt_fit$posterior[,1] < 0.5) + 1

# To be compared with its "extreme" version: k-means (with 2 centers)
?kmeans
km_fit <- kmeans( faithful$waiting, centers = c(45, 55) )
names(km_fit)
# Compute the proportion, center and standard deviation for each cluster,
list(p     = km_fit$size/sum(km_fit$size), 
     mu    = as.vector(km_fit$centers), 
     sigma = sqrt(km_fit$withinss/km_fit$size)
)
#...compare...
list(p = mt_fit$lambda, mu = mt_fit$mu, sigma = mt_fit$sigma)

# Model Based :: MoG ------------------------------------------------------

## Geyser Data ---
suppressMessages(require(plotly, quietly = T))
?MASS::geyser
dd <- MASS::geyser
#kde2d to make the kernel density estimation in 2D
kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
p <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
p

# Mixture of Gaussian
suppressMessages(require(mclust, quietly = T))
?mclust

mod1 = Mclust(dd)
mod1
summary(mod1)
mod1$classification

# Default plot
class(mod1)
?plot.Mclust

par(mfrow = c(2,2))
plot(mod1, what = "uncertainty")
plot(mod1, what = "classification")
plot(mod1, what = "density")
plot(mod1, what = "BIC")
par(mfrow = c(1,1))

## Circle ---

# Noiseless
set.seed(123)
n      <- 500
radius <- 2
angle  <- runif(n)*pi*2
x      <- cos(angle)*radius
y      <- sin(angle)*radius

# Jitter around
SD <- .15
xn <- x + rnorm(n, sd = SD)
yn <- y + rnorm(n, sd = SD)
XX <- cbind(xn, yn)
plot(x,y, pch = 19, col = "black", asp = 1, cex = .5)
points(xn, yn, pch = 19, col = "red", cex = .5)

mod2 = Mclust(cbind(xn,yn))
mod2
summary(mod2)

par(mfrow = c(2,2))
plot(mod2, what = "uncertainty")
plot(mod2, what = "classification")
plot(mod2, what = "density")
plot(mod2, what = "BIC")
par(mfrow = c(1,1))

