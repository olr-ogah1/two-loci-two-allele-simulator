haplotype_simulate <- function(initial.freq=c(0.25, 0.25, 0.25, 0.25), N=pop_1, c=0.5, t=100, replicates = 1000){
  #N is a vector of the population size from the past to the most recent
  #FIRST CHECK ALL INPUTS ARE VALID
  problems <- 0
  if (replicates != round(replicates) || replicates == 0){
    print("Check replicates (value has to be a non-zero integer)")
    problems <- 1 + problems
  }
  if (t!= round(t) || t ==0 ) {
    print("Check t (value has to be a non-zero integer)")
    problems <- 1 + problems
  }
  if (length(N) != t) {
    print("Check N and t (there has to be an equivalent element of N for every time point)")
    problems <- 1 + problems
  }
  if (sum(initial.freq) != 1){
    print("Check initial frequencies (all 4 frequencies have to equal 1)")
    problems <- 1 + problems
  }
  if (c > 0.5 || c < 0) {
    print("Check c (c needs to be between 0 and 0.5)")
    problems <- 1 + problems
  }
  if (all(N != round(N)) || any(N ==0)) {
    print("Check vector N (all values need to be integers)")
    problems <- 1 + problems
  }
  if (problems == 0){
    print("There are 0 errors. The simulator is running")
    rsquared_total <- matrix(NA, nrow = t , ncol = replicates) #Create a matrix to store all the rsqaured values 
    replicate_number = 1
    while (replicate_number <= replicates) {
      # CREATE SOME (EMPTY) VECTORS TO STORE OUR OUTPUTS
      pAB<-rep(NA, t+1)
      pAb<-rep(NA, t+1)
      paB<-rep(NA, t+1)
      pab<-rep(NA, t+1)
      D<-rep(NA, t+1)
      rsquared <- rep(NA, t)
      # FILL THESE VECTORS WITH INITIAL VALUES
      pAB[1]<-initial.freq[1]
      pAb[1]<-initial.freq[2]
      paB[1]<-initial.freq[3]
      pab[1]<-initial.freq[4]
      D[1]<-pAB[1]*pab[1]-pAb[1]*paB[1]
      rsquared[1] <- (D[1]/(sqrt((pAB[1] + pAb[1]) * (paB[1] + pab[1]) * (pAB[1]+paB[1]) * (pAb[1]+pab[1]))))^2
      # FOR LOOP FOR EACH GENERATION - PROPAGATION. NEED TO CONSIDER BOTH DRIFT AND RECOMBINATION
      for (i in 1:t) {
        # CALCULATE THE GAMETIC FREQ IN THE GAMETE POOL
        # WILL BE ERASED AND RE-WRITTEN IN EACH TIME STEP
        gameteAB<-pAB[i]-c*D[i]
        gameteAb<-pAb[i]+c*D[i]
        gameteaB<-paB[i]+c*D[i]
        gameteab<-pab[i]-c*D[i]
        # DRIFT, RANDOM SAMPLING FROM THE GAMETE POOL
        individuals_in_gen<-rmultinom(1, size=2*N[i], prob=c(gameteAB, gameteAb, gameteaB, gameteab))
        pAB[i+1]<-individuals_in_gen[1]/(2*N[i])
        pAb[i+1]<-individuals_in_gen[2]/(2*N[i])
        paB[i+1]<-individuals_in_gen[3]/(2*N[i])
        pab[i+1]<-individuals_in_gen[4]/(2*N[i])
        if (pAB[i+1] <= 0.05 || pAb[i+1] <= 0.05 || paB[i+1]<= 0.05 || pab[i+1]<= 0.05 || pAB[i+1]>= 0.95 || pAb[i+1]>= 0.95 || paB[i+1]>= 0.95 || pab[i+1]>= 0.95){
          rsquared <- 2
          break
        } else {
          #CALCULATE R^2 FOR THIS GENERATION
          rsquared[i] <- (D[i]/(sqrt((pAB[i] + pAb[i])* (paB[i] + pab[i]) * (pAB[i]+paB[i]) * (pAb[i]+pab[i]))))^2
          # CALCULATE THE D IN THE NEXT GENERATION
          D[i+1]<-pAB[i+1]*pab[i+1]-pAb[i+1]*paB[i+1]
        }
      }
      #FILL IN THE MATRIX WITH R-SQUARED VALUES FOR EACH REPLICATE SEQUENTIALLY 
      if (rsquared[1] != 2){
        rsquared_total[,replicate_number] <- rsquared
        replicate_number = replicate_number + 1
      }
    }
    return(rsquared_total) 
  } 
}
