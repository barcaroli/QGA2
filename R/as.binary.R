#---------------------------
# FROM DECIMAL TO BINARY                    
#---------------------------  
as.binary <- function(number,n) {
  bin <- rep(NA,n)
  i = n
  for (i in c(n:1)) {
    digit <- number %% 2
    number <- floor(number / 2)
    bin[i] <- digit
  }
  return(bin)
}
