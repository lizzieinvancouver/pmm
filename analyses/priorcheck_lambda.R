## Code for visualizing the prior distribution of Pagel's lambda
## Can be used to set reasonable values for alpha and beta in the lambda priors

## Set alpha, beta
alpha <- 14
beta <- 6
## Range of values to plot over from 0 to 1
xrange <- seq(0, 1, length.out = 100)
## Make plot
plot(dbeta(x = xrange, shape1 = alpha, shape2 = beta) ~ xrange, type = "l")
