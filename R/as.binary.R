#---------------------------
# FROM DECIMAL TO BINARY                    
#---------------------------  
as.binary <- function(number, n) {
  # Use integer division and modular arithmetic to extract bits
  number <- as.integer(number)
  bits <- rev(as.integer(intToBits(number))[1:n])
  return(bits)
}
