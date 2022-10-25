#------------------
# REPAIR PROCEDURE                     
#------------------
repair <- function(chromosome) {
  diff = 2^Genome_el - nstrat
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(Genome_el,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:Genome_el)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
      }
    }
    solution <- solution + 1
    table(solution)
    sum(table(solution))
    t <- as.numeric(table(solution))
    length(t)
    if (length(t) > nstrat) { 
      solution[solution %in% c(which(t==min(t))[1]:length(t)) & !(solution %in% c(1:diff))] <- solution[solution %in% c(which(t==min(t))[1]:length(t)) & !(solution %in% c(1:diff))] - diff
    }
    a = array(c(1:genomeLength),c(Genome_el,Genome))
    for (x in c(1:Genome)) {
      y1 = a[1,x]
      y2 = a[Genome_el,x]
      chromosome[i,c(y1:y2)] <- as.binary(solution[x]-1,n=Genome_el)
    }
  }  
  return(chromosome)
}

