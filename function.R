n <- 128
x <- c(1:n)
y <- rep(NA,n)
for (i in c(1:n)) {
  y[i] <- (x[i]-n/2)^2
  y[i] <- y[i] + runif(1,min=-2*n,max=2*n)
}
plot(x,y)
which(y[x]==min(y))
