####################
### CODE ARTICLE ###
####################

################################################################################

####################
### 0. Libraries ###
####################

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
set.seed(1234)

####################
### A) Functions ###
####################

##### 1. Matrix generation #####

#Squared matrix simulation for individual sequence
#Rectangular matrix simulation for dyadic sequences
#Simulation using uniform distribution
#Dyadic sequences: 3 patterns (actor-partner, actor only, partner only)

indMatrix <- function(states){
  TP <- matrix(NA, nrow=states, ncol=states)
  for(i in 1:nrow(TP)){
    uniform <- runif(states)
    prob <- (uniform/sum(uniform))
    TP[i,] <- prob
  }
  if(!all(near(rowSums(TP),1))){
    stop("Error in the row of the probability transition matrix")
  }
  if(any(TP==0)){
    stop("Presence of 0 in the probability transition matrix")
  }
  return(TP)
}

apimMatrix <- function(states, type = c("APM", "AM", "PM")){
  TP <- matrix(NA, nrow=states*states, ncol=states)
  if(type == "APM"){
    for(i in 1:nrow(TP)){
      uniform <- runif(states)
      prob <- (uniform/sum(uniform))
      TP[i,] <- prob
    }
  }
  if(type == "AM"){
    for(i in 1:states){
      uniform <- runif(states)
      prob <- (uniform/(sum(uniform)))
      for(j in 1:states){
        TP[(j+states*(i-1)),] <- prob
      }
    }
  }
  if(type == "PM"){
    for(i in 1:states){
      uniform <- runif(states)
      prob <- (uniform/sum(uniform))
      for(j in 1:states){
        TP[(i+states*(j-1)),] <- prob
      }
    }
  } 
  if(!all(near(rowSums(TP),1))){
    stop("Error in the row of the probability transition matrix")
  }
  if(any(TP==0)){
    stop("Presence of 0 in the probability transition matrix")
  }
  return(TP)
}



##### 2. Find a row #####

#Preliminary function for the generation of the sequence for dyadic case
#Given state for one member, find the corresponding row in the partner matrix

findRow <- function(stateFM, stateSM, nbcol){ 
  row <- (((stateFM-1)*nbcol)+stateSM)
  return(row)
}



##### 3. Sequence generation #####

#Two approaches for the sequences generation: with or without matrix in input
#"Alternative" does not simulate the matrix -> must be given 
#"Sim" simulates the matrix and generates the sequence 

### Discrete MC simulation

indSimAlt <- function(measure, matSim, statesInitial = NA){
  TP <- matSim
  if(!near(ncol(TP), nrow(TP))){
    stop("Error in the dimension of the probability transition matrix")
  }
  chain <- rep(0, measure)
  if(!any(is.na(statesInitial))){
    chain[1] <- statesInitial
  }
  else {
    row <- sample(1:dim(TP)[1], 1)
    rowProb <- TP[row,]
    initialMultinomial <- rmultinom(1,1,rowProb)
    chain[1] <- which(initialMultinomial==1)
  }
  for(t in 2:measure){
    probState <- TP[chain[t-1],]
    multinomial <- rmultinom(1,1,probState)
    chain[t] <- which(multinomial==1)
  }
  return(chain)
}

indSim <- function(states, measure, statesInitial = NA){
  TP <- indMatrix(states)
  if(!near(ncol(TP), nrow(TP))){
    stop("Error in the dimension of the probability transition matrix")
  }
  chain <- rep(0, measure)
  if(!any(is.na(statesInitial))){
    chain[1] <- statesInitial
  }
  else {
    row <- sample(1:dim(TP)[1], 1)
    rowProb <- TP[row,]
    initialMultinomial <- rmultinom(1,1,rowProb)
    chain[1] <- which(initialMultinomial==1)
  }
  for(t in 2:measure){
    probState <- TP[chain[t-1],]
    multinomial <- rmultinom(1,1,probState)
    chain[t] <- which(multinomial==1)
  }
  return(list(TP, chain))
}

apimSimAlt <- function(measure, matSimFM, matSimSM, statesInitial = NA){
  TPX <- matSimFM
  TPY <- matSimSM
  if(!identical(dim(TPX),dim(TPY))){
    stop("Error in the dimension of the probability transition matrices")
  }
  chainX <- rep(0, measure) 
  chainY <- rep(0, measure)
  if(!any(is.na(statesInitial))){
    chainX[1] <- statesInitial[1]
    chainY[1] <- statesInitial[2]
  }
  else {
    initialX <- sample(1:ncol(TPX), 1)
    initialY <- sample(1:ncol(TPY), 1)
    rowXProb <- TPX[findRow(initialX, initialY, ncol(TPX)),]
    rowYProb <- TPY[findRow(initialY, initialX, ncol(TPY)),]
    initialXMultinomial <- rmultinom(1,1,rowXProb)
    initialYMultinomial <- rmultinom(1,1,rowYProb)
    chainX[1] <- which(initialXMultinomial==1)
    chainY[1] <- which(initialYMultinomial==1)
  }
  for(t in 2:measure){
    previousX <- chainX[t-1]
    previousY <- chainY[t-1]
    probXState <- TPX[findRow(previousX, previousY, ncol(TPX)),] 
    probYState <- TPY[findRow(previousY, previousX, ncol(TPY)),] 
    multinomialX <- rmultinom(1,1,probXState)
    multinomialY <- rmultinom(1,1,probYState) 
    chainX[t] <- which(multinomialX==1) 
    chainY[t] <- which(multinomialY==1) 
  }
  return(list(chainX, chainY))
}

apimSim <- function(states, measure, statesInitial = NA, type = c("APM", "AM", "PM")){
  TPX <- apimMatrix(states = states, type = type)
  TPY <- apimMatrix(states = states, type = type)
  if(!identical(dim(TPX),dim(TPY))){
    stop("Error in the dimension of the probability transition matrices")
  }
  chainX <- rep(0, measure) 
  chainY <- rep(0, measure)
  if(!any(is.na(statesInitial))){
    chainX[1] <- statesInitial[1]
    chainY[1] <- statesInitial[2]
  }
  else {
    initialX <- sample(1:ncol(TPX), 1)
    initialY <- sample(1:ncol(TPY), 1)
    rowXProb <- TPX[findRow(initialX, initialY, ncol(TPX)),]
    rowYProb <- TPY[findRow(initialY, initialX, ncol(TPY)),]
    initialXMultinomial <- rmultinom(1,1,rowXProb)
    initialYMultinomial <- rmultinom(1,1,rowYProb)
    chainX[1] <- which(initialXMultinomial==1)
    chainY[1] <- which(initialYMultinomial==1)
  }
  for(t in 2:measure){
    previousX <- chainX[t-1]
    previousY <- chainY[t-1]
    probXState <- TPX[findRow(previousX, previousY, ncol(TPX)),] 
    probYState <- TPY[findRow(previousX, previousY, ncol(TPY)),] 
    multinomialX <- rmultinom(1,1,probXState)
    multinomialY <- rmultinom(1,1,probYState) 
    chainX[t] <- which(multinomialX==1) 
    chainY[t] <- which(multinomialY==1) 
  }
  return(list(TPX, TPY, chainX, chainY))
}



##### 4. Estimation #####

#To estimate the transition probabilities matrix, need for the counts 
#If there is no count, we replace "Nan" with 0 -> no transition
#Estimation of the transition probabilities matrix using maximum likelihood 

countEmp <- function(states, chainFM, chainSM = NULL){
  chainCount <- (length(chainFM)-1)
  if(is.null(chainSM)){
    count <- matrix(0, nrow = states, ncol = states)
    for(i in 1:chainCount){
      column <- chainFM[i+1]
      row <- chainFM[i]
      count[row, column] <- (count[row, column]+1)
    }
  } else {
    count <- matrix(0, nrow = states*states, ncol = states)
    for(i in 1:chainCount){
      column <- chainFM[i+1]
      behaviorFM <- chainFM[i]
      behaviorSM <- chainSM[i]
      row <- ((1+states*(behaviorFM-1))+(behaviorSM-1))
      count[row, column] <- (count[row, column]+1)
    }
  }
  if(sum(count)!=(length(chainFM)-1)){
    stop("Error in the count")
  }
  return(count) 
}

applyZero <- function(matrix){
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(is.nan(matrix[i,j])){
        matrix[i,j] <- 0
      }
    }
  }
  return(matrix)
}

mleEstimation <- function(countMat){
  estimate <- matrix(0, nrow = nrow(countMat), ncol = ncol(countMat))
  rowSum <- rowSums(countMat)
  for(i in 1:nrow(countMat)){
    vectorProb <- (countMat[i,]/rowSum[i])
    estimate[i,] <- vectorProb
  }
  if(any(is.nan(estimate))){
    estimate <- applyZero(estimate)
  }
  return(estimate)
}



##### 5. Hypothesis testing #####

#All the test are done using chi square distance
#Need of the empirical and theoretical count matrix 
#The theoretical count matrix can be computed according to restrictions or not
#The chi squared test of GOF can be done for each sequence or for both
#The LRT test can be done for each sequence or for both 
#GOF = squaredTest, LRT = difftest
#1 sequence = local, 2 sequences = global

chisquaredDist <- function(population, empirical){
  chidist <- 0
  for(i in 1:nrow(population)){
    for(j in 1:ncol(population)){
      numerator <- ((empirical[i,j]-population[i,j])^2)
      denominator <- population[i,j]
      if(denominator == 0){
        ratio <- 0
      } else{
        ratio <- (numerator/denominator)
      }
      chidist <- (chidist+ratio)
    }
  }
  return(chidist)
}

countTheo <- function(lambda, countEmp, condition = c("U", "C"), type = c("AM","PM")){
  gamma <- rowSums(countEmp)
  states <- ncol(countEmp)
  xi <- vector()
  int <- matrix(NA, ncol = ncol(countEmp), nrow = nrow(countEmp))
  neta <- matrix(NA, nrow=ncol(countEmp), ncol=ncol(countEmp))
  countTheo <- matrix(NA, ncol = ncol(countEmp), nrow = nrow(countEmp))
  if(condition == "U"){
    for(i in 1:nrow(countEmp)){
      countTheo[i,] <- (gamma[i]*lambda[i,])
    }
  } else {
    if(type == "AM"){
      sumGroup <- (diag(rep(1, ncol(countEmp))) %x% t(rep(1, ncol(countEmp))) %*% countEmp)
      for(i in 1:states){
        xiValue <- sum(gamma[(1+(i-1)*states):(i*states)])
        xi <- append(xi, xiValue)
      }
      for(i in 1:length(xi)){
        neta[i,] <- sumGroup[i,]/xi[i]
      }
      for(i in 1:ncol(neta)){
        for(j in 1:ncol(neta)){
          int[(j+states*(i-1)),] <- neta[i,]
        }
      }
      for(i in 1:nrow(countEmp)){
        countTheo[i,] <- (int[i,]*gamma[i])
      }
    } else {
      sumGroup <- ((do.call(cbind, replicate(states, diag(rep(1, ncol(countEmp))), simplify = FALSE)))  %*% countEmp)
      for(i in 1:states){
        num <- vector()
        for(j in 1:states){
          pos <- gamma[i + states*(j-1)]
          num <- append(num, pos)
        }
        xiValue <- sum(num)
        xi <- append(xi, xiValue)
      }
      for(i in 1:length(xi)){
        neta[i,] <- sumGroup[i,]/xi[i]
      }
      int <- do.call(rbind, replicate(states, neta, simplify=FALSE)) 
      for(i in 1:nrow(countEmp)){
        countTheo[i,] <- (int[i,]*gamma[i])
      }
    }
  }
  return(countTheo)
}

accuracyLocal <- function(population, empirical){
  method      <- "Chisquared test"
  dataName   <- "Observed vs Estimated"
  alternative <- "The model does not fit the data"
  khi2 <- chisquaredDist(population = population, empirical = empirical)
  degree <- ((ncol(population)-1)*nrow(population))
  pValue <- pchisq(q=khi2, df=degree, lower.tail = F)
  names(khi2) <- "X-squared"
  names(degree) <- "df"
  TEST        <- list(method = method, data.name = dataName,
                      parameter = degree, alternative = alternative, statistic = khi2, p.value = pValue)
  class(TEST) <- "htest"
  TEST
}

accuracyGlobal <- function(populationFM, empiricalFM, populationSM, empiricalSM){
  method      <- "Chisquared test"
  dataName   <- "Observed vs Estimated at Level 1 "
  alternative <- "The model does not fit the data"
  khiFM <- chisquaredDist(population = populationFM, empirical = empiricalFM)
  khiSM <- chisquaredDist(population = populationSM, empirical = empiricalSM)
  khi2 <- (khiFM + khiSM)
  degree <- (2*((ncol(populationFM)-1)*nrow(populationFM)))
  pValue <- pchisq(q=khi2, df=degree, lower.tail = F)
  names(khi2) <- "X-squared"
  names(degree) <- "df"
  TEST        <- list(method = method, data.name = dataName,
                      parameter = degree, alternative = alternative, statistic = khi2, p.value = pValue)
  class(TEST) <- "htest"
  TEST
}

lrtLocal <- function(population, empirical){
  method      <- "Chisquared test"
  dataName   <- "Observed vs Estimated"
  alternative <- "The model does not fit the data"
  khi2 <- chisquaredDist(population = population, empirical = empirical)
  degree <- (ncol(population)*(ncol(population)-1)^2)
  pValue <- pchisq(q=khi2, df=degree, lower.tail = F)
  names(khi2) <- "X-squared"
  names(degree) <- "df"
  TEST        <- list(method = method, data.name = dataName,
                      parameter = degree, alternative = alternative, statistic = khi2, p.value = pValue)
  class(TEST) <- "htest"
  TEST
}

lrtGlobal <- function(populationFM, empiricalFM, populationSM, empiricalSM){
  method      <- "Chisquared test"
  dataName   <- "Observed vs Estimated"
  alternative <- "The model does not fit the data"
  khiFM <- chisquaredDist(population = populationFM, empirical = empiricalFM)
  khiSM <- chisquaredDist(population = populationSM, empirical = empiricalSM)
  khi2 <- (khiFM + khiSM)
  degree <- (2*(ncol(populationFM)*(ncol(populationFM)-1)^2))
  pValue <- pchisq(q=khi2, df=degree, lower.tail = F)
  names(khi2) <- "X-squared"
  names(degree) <- "df"
  TEST        <- list(method = method, data.name = dataName,
                      parameter = degree, alternative = alternative, statistic = khi2, p.value = pValue)
  class(TEST) <- "htest"
  TEST
}



##### 6. Accuracy #####

#The accuracy measure is the root mean squared error 
#Compute the RMSE according to simulated sequences 
#Compute the RMSE and compare to given quantile

rmse <- function(lambdaTheo, lambdaEst, type = c("ind", "dyad")){
  s <- ncol(lambdaTheo)
  if(type == "ind"){
    denominator <- (s^2)
    it <- s
  } else {
    denominator <- (s^3)
    it <- (s^2)
  }
  numerator <- 0
  for(i in 1:it){
    for(j in 1:s){
      numerator <- (numerator + (sum((lambdaEst[i,j]-lambdaTheo[i,j])^2)))
    }
  }
  accuracy <- (sqrt(numerator/denominator))
  return(accuracy)
}

errorAccuracy <- function(nbsim, states, lgth, model, person = c("ind", "dyad")){
  TAB <- NULL
  for(n in 1:nbsim){
    lambdaSimX <- apimMatrix(states = states, type = model) #matrix wrt model 
    lambdaSimY <- apimMatrix(states = states, type = model) #matrix wrt model 
    ERROR <- NULL
    for(k in 1:length(lgth)){
      simulation <- apimSimAlt(measure = lgth[k], matSimFM = lambdaSimX, matSimSM = lambdaSimY) #simulation
      chainX <- simulation[[1]] #chain for the first member (actor)
      chainY <- simulation[[2]] #chain for the second member (partner)
      countEmpX <- countEmp(states = states, chainFM =  chainX, chainSM = chainY) #empirical count matrix 
      countEmpY <- countEmp(states = states, chainFM =  chainY, chainSM = chainX) #empirical count matrix 
      lambdaEstX <- mleEstimation(countMat = countEmpX)
      lambdaEstY <- mleEstimation(countMat = countEmpY)
      errorFM <- rmse(lambdaSimX, lambdaEstX, type = person) #error 
      errorSM <- rmse(lambdaSimY, lambdaEstY, type = person) #error
      error <- sqrt((errorFM^2 + errorSM^2)/2)
      ERROR <- c(ERROR, errorFM, errorSM)
    }
    TAB <- rbind(TAB, ERROR)
  }
  return(TAB)
}   

testAccuracy<-function(accuFM,accuSM, threshold, quantile){
  nbsim<-length(accuFM)
  count<-0
  #count element
  for(i in 1:nbsim){
    if(accuFM[i]<threshold && accuSM[i]<threshold){
      count<-count+1
    }
  }
  if(count>(quantile*nbsim)){
    return(TRUE)    
  }else{
    return(FALSE)
  }
}





################################
### B) Length and Simulation ###
################################

##### 1. Parameters #####

states <- c(2,3) #2 states and 3 states sequences
lgth <- c(20, 50, 100, 200, 500, 1000, 2000) #number of points 
nbsim <-c(100, 500, 1000, 10000) #number of simulations

#names for the database
twothreeNames = c("two", "three")
dataNames<-vector()
for(s in 1:length(states)){
  for(n in 1:length(nbsim)){
    dataNames<-append(dataNames,paste(twothreeNames[s],toString(nbsim[n]),sep = ""))
  }
}



##### 2. Length analysis #####

##### 2.1.Database creation #####

#Database with error given states, number of points and number of simulation 
#Variables: error (RMSE), length, member of the dyad (FM, SM)
#Creation of a interaction variable between the length and the member (plots)

for(s in 1:length(states)){
  for(n in 1:length(nbsim)){
    data <- NULL
    #print(dataNames[n+(s-1)*4])
    res <- errorAccuracy(nbsim[n], states[s], lgth, "APM", "dyad")
    data <- data.frame(
      Error = c(res),
      Length = rep(lgth, rep(2*nbsim[n], length(lgth))),
      Member = rep(rep(c("1", "2"), length(lgth)), rep(nbsim[n], 2*length(lgth)))
    )
    data$int <- interaction(data$Length, data$Member)
    assign(dataNames[n+(s-1)*4],data)
  }
}

##### 2.2 Plots analysis #####

plotTwo100 <- ggplot(two100, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "100 simulations")

plotTwo500 <- ggplot(two500, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "500 simulations")

plotTwo1000 <- ggplot(two1000, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "1000 simulations")

plotTwo10000 <- ggplot(two10000, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0,0.8) +
  scale_x_discrete(labels = lgth) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "10000 simulations")


mainBox <- textGrob("RMSE boxplot for 2 states",gp=gpar(fontsize=20,font=3))
grid.arrange(plotTwo100,plotTwo500,plotTwo1000,plotTwo10000, ncol=2, nrow=2, top = mainBox)
#dev.print(pdf, file = "RMSE2box.pdf") 
dev.print(png, file = "RMSE2box.png", width = 1000, height = 800) 


plotThree100 <- ggplot(three100, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "100 simulations")

plotThree500 <- ggplot(three500, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "500 simulations")

plotThree1000 <- ggplot(three1000, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "1000 simulations")

plotThree10000 <- ggplot(three10000, aes(x=int, y=Error, fill=Member)) + 
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  scale_x_discrete(labels = lgth) +
  ylim(0,0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Length",
       title = "10000 simulations")


mainBox <- textGrob("RMSE boxplot for 3 states",gp=gpar(fontsize=20,font=3))
grid.arrange(plotThree100,plotThree500,plotThree1000,plotThree10000, ncol=2, nrow=2, top = mainBox)
#dev.print(pdf, file = "RMSE3box.pdf", width = 10, height = 7) 
dev.print(png, file = "RMSE3box.png", width = 1000, height = 800) 


##### 3. Dyadic simulations #####
##### 3.1 Parameters #####

lgth <- 1000 #fix number of points 
dyad <- 500 #number of dyads

#names for the database and the variables 
dyadId <- vector()
for(i in 1:dyad){
  each <- rep(i, (2*lgth))
  dyadId <- append(dyadId, each)   
}
memberId <- rep(c(rep("1", lgth),rep("2", lgth)), dyad)
measurement <- rep(c(1:lgth),2*dyad)

##### 3.2 General simulations #####

#Simulation of two databases: APM, PATTERN
#Database with only Actor-Partner pattern 
#Database with Actor-Partner, Actor-only, Partner-only patterns 
#-> looks like an empirical case 

### A. APM simulation 

dataNames<-vector()
for(s in 1:length(states)){
  dataNames <- append(dataNames, paste(twothreeNames[s],"dyadAPM", sep = ""))
}

for(s in 1:length(states)){
  sequence <- vector()
  for(i in 1:dyad){
    simulation <- apimSim(states = states[s], measure = lgth, type = "APM")
    chainFM <- simulation[[3]]
    chainSM <- simulation[[4]]
    chain <- c(chainFM, chainSM)
    sequence <- append(sequence, chain)
  }
  data <- data.frame(cbind(dyadId, memberId, measurement, sequence))
  assign(dataNames[s],data)
}

### B. Patterns simulation 

models <- c("APM", "AM", "PM")
dataNames<-vector()
for(s in 1:length(states)){
  dataNames <- append(dataNames, paste(twothreeNames[s],"dyadPATTERN", sep = ""))
}

for(s in 1:length(states)){
  sequence <- vector()
  pattern <- vector()
  for(i in 1:dyad){
    model <- sample(models, 1)
    pattern <- append(pattern, rep(model, (2*lgth)))
    simulation <- apimSim(states = s, measure = lgth, type = model)
    chainFM <- simulation[[3]]
    chainSM <- simulation[[4]]
    chain <- c(chainFM, chainSM)
    sequence <- append(sequence, chain)
  }
  data <- data.frame(cbind(dyadId, memberId, measurement, pattern, sequence))
  assign(dataNames[s],data)
}


##### 3.4 Prototypic simulations #####

#Simulation with prototypic matrices 
#Equiprobable transition matrix 
#Extreme transition matrix 

### A. Equiprobable matrix 

dataNames<-vector()
for(s in 1:length(states)){
  dataNames <- append(dataNames, paste(twothreeNames[s],"dyadEQU", sep = ""))
}

for(s in 1:length(states)){
  lambdaFM <- matrix(rep(1/states[s], states[s]^3), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  lambdaSM <- matrix(rep(1/states[s], states[s]^3), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  sequence <- vector()
  for(i in 1:dyad){
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM)
    chainFM <- simulation[[1]]
    chainSM <- simulation[[2]]
    chain <- c(chainFM, chainSM)
    sequence <- append(sequence, chain)
  }
  data <- data.frame(cbind(dyadId, memberId, measurement, sequence))
  assign(dataNames[s],data)
}

### B. Extreme matrix 

dataNames<-vector()
for(s in 1:length(states)){
  dataNames <- append(dataNames, paste(twothreeNames[s],"dyadEXT", sep = ""))
}

for(s in 1:length(states)){
  lambdaFM <- matrix(c((1-(sum(rep(0.1, (s-1))))), (rep(0.1, (s-1)))), ncol = s, nrow = s^2, byrow = TRUE)
  lambdaSM <- matrix(c((1-(sum(rep(0.1, (s-1))))), (rep(0.1, (s-1)))), ncol = s, nrow = s^2, byrow = TRUE)
  sequence <- vector()
  for(i in 1:dyad){
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM)
    chainFM <- simulation[[1]]
    chainSM <- simulation[[2]]
    chain <- c(chainFM, chainSM)
    sequence <- append(sequence, chain)
  }
  data <- data.frame(cbind(dyadId, memberId, measurement, sequence))
  assign(dataNames[s],data)
}

##############
### C) GOF ###
##############

##### 1. Parameters #####

nbsim <- 1000 #fix number of simulations
quant <- c(0.8, 0.9, 0.95) #quantile for chi squared distribution

members <- c("FM", "SM", "Global") #local and global test 

simNames <- c(1:nbsim)
memNames <- members
statesNames <- c("2states","3states") #states names for table
quantNames <- quant

##### 2. Actor-Partner matrix pattern #####

model <- c("APM")
modNames <- model 
tableNames <- paste("tableGOF", modNames, sep = "")
testNames <- paste("testGOF", modNames, sep = "")

tableau <- array(NA, dim = c(nbsim, length(members), length(states)), 
                 dimnames = list(simNames, memNames, statesNames))
test <- array(NA, dim = c(length(quant), length(members), length(states)), 
              dimnames = list(quantNames, memNames, statesNames))

for(s in 1:length(states)){
  lambdaFM <- apimMatrix(states = states[s], type = model) #matrix wrt model 
  lambdaSM <- apimMatrix(states = states[s], type = model) #matrix wrt model
  for(i in 1:nbsim){ #run the simulation
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM) #simulation
    chainX <- simulation[[1]] #chain for the first member (actor)
    chainY <- simulation[[2]] #chain for the second member (partner)
    countEmpX <- countEmp(states = states[s], chainFM =  chainX, chainSM = chainY) #empirical count matrix 
    countEmpY <- countEmp(states = states[s], chainFM =  chainY, chainSM = chainX) #empirical count matrix 
    countTheoX <- countTheo(lambda = lambdaFM, countEmp = countEmpX, condition = "U")
    countTheoY <- countTheo(lambda = lambdaSM, countEmp = countEmpY, condition = "U")
    compFM <- accuracyLocal(population = countTheoX, empirical = countEmpX) #test 
    compSM <- accuracyLocal(population = countTheoY, empirical = countEmpY) #test
    compGlobal <- accuracyGlobal(populationFM = countTheoX, populationSM = countTheoY, empiricalFM = countEmpX, countEmpY)
    distanceFM <- compFM$statistic[["X-squared"]] #distance
    distanceSM <- compSM$statistic[["X-squared"]] #distance
    distanceGlobal <- compGlobal$statistic[["X-squared"]]
    tableau[i,1,s] <- distanceFM
    tableau[i,2,s] <- distanceSM
    tableau[i,3,s] <- distanceGlobal
  }
}

for(s in 1:length(states)){
  df <- ((states[s])^2)*(states[s]-1)
  for(q in 1:length(quant)){
    test[q,1,s] <- ((sum(tableau[,1,s] > qchisq(quant[q], df)))/nbsim)
    test[q,2,s] <- ((sum(tableau[,2,s] > qchisq(quant[q], df)))/nbsim)
    test[q,3,s] <- ((sum(tableau[,3,s] > qchisq(quant[q], (2*df))))/nbsim)
  }
}

assign(tableNames, tableau)
assign(testNames, test)

##### 3. Prototypic patterns #####

##### 3.1 Equiprobable matrix #####

tableNames <- paste("tableGOF", "EQU", sep = "")
testNames <- paste("testGOF", "EQU", sep = "")

tableau <- array(NA, dim = c(nbsim, length(members), length(states)), 
                 dimnames = list(simNames, memNames, statesNames))
test <- array(NA, dim = c(length(quant), length(members), length(states)), 
              dimnames = list(quantNames, memNames, statesNames))

for(s in 1:length(states)){
  lambdaFM <- matrix(rep(1/states[s], states[s]^3), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  lambdaSM <- matrix(rep(1/states[s], states[s]^3), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  for(i in 1:nbsim){ 
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM) 
    chainX <- simulation[[1]] 
    chainY <- simulation[[2]] 
    countEmpX <- countEmp(states = states[s], chainFM =  chainX, chainSM = chainY) 
    countEmpY <- countEmp(states = states[s], chainFM =  chainY, chainSM = chainX) 
    countTheoX <- countTheo(lambda = lambdaFM, countEmp = countEmpX, condition = "U")
    countTheoY <- countTheo(lambda = lambdaSM, countEmp = countEmpY, condition = "U")
    compFM <- accuracyLocal(population = countTheoX, empirical = countEmpX) 
    compSM <- accuracyLocal(population = countTheoY, empirical = countEmpY) 
    compGlobal <- accuracyGlobal(populationFM = countTheoX, populationSM = countTheoY, empiricalFM = countEmpX, countEmpY)
    distanceFM <- compFM$statistic[["X-squared"]]
    distanceSM <- compSM$statistic[["X-squared"]] 
    distanceGlobal <- compGlobal$statistic[["X-squared"]]
    tableau[i,1,s] <- distanceFM
    tableau[i,2,s] <- distanceSM
    tableau[i,3,s] <- distanceGlobal
  }
}

for(s in 1:length(states)){
  df <- ((states[s])^2)*(states[s]-1)
  for(q in 1:length(quant)){
    test[q,1,s] <- ((sum(tableau[,1,s] > qchisq(quant[q], df)))/nbsim)
    test[q,2,s] <- ((sum(tableau[,2,s] > qchisq(quant[q], df)))/nbsim)
    test[q,3,s] <- ((sum(tableau[,3,s] > qchisq(quant[q], (2*df))))/nbsim)
  }
}

assign(tableNames, tableau)
assign(testNames, test)

##### 3.2 Extreme matrix #####

tableNames <- paste("tableGOF", "EXT", sep = "")
testNames <- paste("testGOF", "EXT", sep = "")

tableau <- array(NA, dim = c(nbsim, length(members), length(states)), 
                 dimnames = list(simNames, memNames, statesNames))
test <- array(NA, dim = c(length(quant), length(members), length(states)), 
              dimnames = list(quantNames, memNames, statesNames))

for(s in 1:length(states)){
  lambdaFM <- matrix(c((1-(sum(rep(0.1, (states[s]-1))))), (rep(0.1, (states[s]-1)))), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  lambdaSM <- matrix(c((1-(sum(rep(0.1, (states[s]-1))))), (rep(0.1, (states[s]-1)))), ncol = states[s], nrow = states[s]^2, byrow = TRUE)
  for(i in 1:nbsim){ 
    simulation <- apiSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM) 
    chainX <- simulation[[1]] 
    chainY <- simulation[[2]] 
    countEmpX <- countEmp(states = states[s], chainFM =  chainX, chainSM = chainY) 
    countEmpY <- countEmp(states = states[s], chainFM =  chainY, chainSM = chainX) 
    countTheoX <- countTheo(lambda = lambdaFM, countEmp = countEmpX, condition = "U")
    countTheoY <- countTheo(lambda = lambdaSM, countEmp = countEmpY, condition = "U")
    compFM <- accuracyLocal(population = countTheoX, empirical = countEmpX) 
    compSM <- accuracyLocal(population = countTheoY, empirical = countEmpY) 
    compGlobal <- accuracyGlobal(populationFM = countTheoX, populationSM = countTheoY, empiricalFM = countEmpX, countEmpY)
    distanceFM <- compFM$statistic[["X-squared"]]
    distanceSM <- compSM$statistic[["X-squared"]] 
    distanceGlobal <- compGlobal$statistic[["X-squared"]]
    tableau[i,1,s] <- distanceFM
    tableau[i,2,s] <- distanceSM
    tableau[i,3,s] <- distanceGlobal
  }
}

for(s in 1:length(states)){
  df <- ((states[s])^2)*(states[s]-1)
  for(q in 1:length(quant)){
    test[q,1,s] <- ((sum(tableau[,1,s] > qchisq(quant[q], df)))/nbsim)
    test[q,2,s] <- ((sum(tableau[,2,s] > qchisq(quant[q], df)))/nbsim)
    test[q,3,s] <- ((sum(tableau[,3,s] > qchisq(quant[q], (2*df))))/nbsim)
  }
}

assign(tableNames, tableau)
assign(testNames, test)



##############
### D) LRT ###
##############

#Parameters are the same as for GOF

##### 1. APM vs AM #####

#Run a LRT between an APM and an AM
#df are the differences between df(APM) and df(AM)

tableNames <- paste("tableLRT", "AM", sep = "")
testNames <- paste("testLRT", "AM", sep = "")

tableau <- array(NA, dim = c(nbsim, length(members), length(states)), 
                 dimnames = list(simNames, memNames, statesNames))
test <- array(NA, dim = c(length(quant), length(members), length(states)), 
              dimnames = list(quantNames, memNames, statesNames))

for(s in 1:length(states)){
  lambdaFM <- apimMatrix(states = states[s], type = "APM")
  lambdaSM <- apimMatrix(states = states[s], type = "APM")
  for(i in 1:nbsim){ 
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM) 
    chainX <- simulation[[1]] 
    chainY <- simulation[[2]] 
    countEmpX <- countEmp(states = states[s], chainFM =  chainX, chainSM = chainY) 
    countEmpY <- countEmp(states = states[s], chainFM =  chainY, chainSM = chainX) 
    countTheoX <- countTheo(lambda = lambdaFM, countEmp = countEmpX, condition = "C", type = "AM")
    countTheoY <- countTheo(lambda = lambdaSM, countEmp = countEmpY, condition = "C", type = "AM")
    compFM <- lrtLocal(population = countTheoX, empirical = countEmpX) 
    compSM <- lrtLocal(population = countTheoY, empirical = countEmpY)
    comp <- lrtGlobal(populationFM = countTheoX, empiricalFM = countEmpX, populationSM = countTheoY, empiricalSM = countEmpY)
    distanceFM <- compFM$statistic[["X-squared"]]
    distanceSM <- compSM$statistic[["X-squared"]] 
    distance <- comp$statistic[["X-squared"]]
    tableau[i,1,s] <- distanceFM
    tableau[i,2,s] <- distanceSM
    tableau[i,3,s] <- distance
  }
}

for(s in 1:length(states)){
  df <- ((states[s]))*((states[s]-1)^2)
  for(q in 1:length(quant)){
    test[q,1,s] <- ((sum(tableau[,1,s] > qchisq(quant[q], df)))/nbsim)
    test[q,2,s] <- ((sum(tableau[,2,s] > qchisq(quant[q], df)))/nbsim)
    test[q,3,s] <- ((sum(tableau[,3,s] > qchisq(quant[q], (2*df))))/nbsim)
  }
}

assign(tableNames, tableau)
assign(testNames, test)

##### 2. APM vs PM #####

#Run a LRT between an APM and an PM
#df are the differences between df(PPM) and df(PM)

tableNames <- paste("tableLRT", "PM", sep = "")
testNames <- paste("testLRT", "PM", sep = "")

tableau <- array(NA, dim = c(nbsim, length(members), length(states)), 
                 dimnames = list(simNames, memNames, statesNames))
test <- array(NA, dim = c(length(quant), length(members), length(states)), 
              dimnames = list(quantNames, memNames, statesNames))

for(s in 1:length(states)){
  lambdaFM <- apimMatrix(states = states[s], type = "APM")
  lambdaSM <- apimMatrix(states = states[s], type = "APM")
  for(i in 1:nbsim){ 
    simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaFM, matSimSM = lambdaSM) 
    chainX <- simulation[[1]] 
    chainY <- simulation[[2]] 
    countEmpX <- countEmp(states = states[s], chainFM =  chainX, chainSM = chainY) 
    countEmpY <- countEmp(states = states[s], chainFM =  chainY, chainSM = chainX) 
    countTheoX <- countTheo(lambda = lambdaFM, countEmp = countEmpX, condition = "C", type = "PM")
    countTheoY <- countTheo(lambda = lambdaSM, countEmp = countEmpY, condition = "C", type = "PM")
    compFM <- lrtLocal(population = countTheoX, empirical = countEmpX) 
    compSM <- lrtLocal(population = countTheoY, empirical = countEmpY)
    comp <- lrtGlobal(populationFM = countTheoX, empiricalFM = countEmpX, populationSM = countTheoY, empiricalSM = countEmpY)
    distanceFM <- compFM$statistic[["X-squared"]]
    distanceSM <- compSM$statistic[["X-squared"]] 
    distance <- comp$statistic[["X-squared"]]
    tableau[i,1,s] <- distanceFM
    tableau[i,2,s] <- distanceSM
    tableau[i,3,s] <- distance
  }
}

for(s in 1:length(states)){
  df <- ((states[s]))*((states[s]-1)^2)
  for(q in 1:length(quant)){
    test[q,1,s] <- ((sum(tableau[,1,s] > qchisq(quant[q], df)))/nbsim)
    test[q,2,s] <- ((sum(tableau[,2,s] > qchisq(quant[q], df)))/nbsim)
    test[q,3,s] <- ((sum(tableau[,3,s] > qchisq(quant[q], (2*df))))/nbsim)
  }
}

assign(tableNames, tableau)
assign(testNames, test)



############################
### E) Empirical Example ###
############################

##### 1. Parameters #####

lgth <- 100
s <- 2

##### 1. General sequences simulation #####

models <- c("APM", "AM", "PM")
members <- c("1", "2")

modNames <- models
memNames <- members

lambda <- array(NA, dim = c(s^2, s, length(models), length(members)),
                dimnames = list(NULL, NULL, modNames, memNames))

sequence <- array(NA, dim = c(lgth, length(members),length(models)), 
                  dimnames = list(NULL, memNames, modNames))

for(j in 1:length(models)){
  simulation <- apimSim(states = s, measure = lgth, statesInitial = NA, type = models[j])
  lambdaX <- simulation[[1]]
  lambdaY <- simulation[[2]]
  chainX <- simulation[[3]]
  chainY <- simulation[[4]]
  lambda[,,j,1]<-lambdaX
  lambda[,,j,2]<-lambdaY
  sequence[,1,j]<-chainX
  sequence[,2,j]<-chainY
}

apm <- data.frame(
  Sequence = as.factor(c(sequence[,1,1],sequence[,2,1])),
  Member = as.factor(c(rep("1",lgth),rep("2",lgth))),
  Model = as.factor(rep("APM", lgth))
)

am <- data.frame(
  Sequence = as.factor(c(sequence[,1,2],sequence[,2,2])),
  Member = as.factor(c(rep("1",lgth),rep("2",lgth))),
  Model = as.factor(rep("AM", lgth))
)

pm <- data.frame(
  Sequence = as.factor(c(sequence[,1,3],sequence[,2,3])),
  Member = as.factor(c(rep("1",lgth),rep("2",lgth))),
  Model = as.factor(rep("PM", lgth))
)

##### 2. Prototypical sequences generation #####

models <- c("EQU", "EXT")
members <- c("1", "2")

modNames <- models
memNames <- members

sequence <- array(NA, dim = c(lgth, length(members),length(models)), 
                  dimnames = list(NULL, memNames, modNames))

lambdaEquX <- matrix(rep(1/s, s^3), ncol = s, nrow = s^2, byrow = TRUE)
lambdaEquY <- matrix(rep(1/s, s^3), ncol = s, nrow = s^2, byrow = TRUE)
simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaEquX, matSimSM = lambdaEquY, statesInitial = NA)
chainX <- simulation[[1]]
chainY <- simulation[[2]]
sequence[,1,1]<-chainX
sequence[,2,1]<-chainY

lambdaExtX <- matrix(c((1-(sum(rep(0.1, (s-1))))), (rep(0.1, (s-1)))), ncol = s, nrow = s^2, byrow = TRUE)
lambdaExtY <- matrix(c((1-(sum(rep(0.1, (s-1))))), (rep(0.1, (s-1)))), ncol = s, nrow = s^2, byrow = TRUE)
simulation <- apimSimAlt(measure = lgth, matSimFM = lambdaExtX, matSimSM = lambdaExtY, statesInitial = NA)
chainX <- simulation[[1]]
chainY <- simulation[[2]]
sequence[,1,2]<-chainX
sequence[,2,2]<-chainY

equ <- data.frame(
  Sequence = as.factor(c(sequence[,1,1],sequence[,2,1])),
  Member = as.factor(c(rep("1",lgth),rep("2",lgth))),
  Model = as.factor(rep("EQU", lgth))
)

ext <- data.frame(
  Sequence = as.factor(c(sequence[,1,2],sequence[,2,2])),
  Member = as.factor(c(rep("1",lgth),rep("2",lgth))),
  Model = as.factor(rep("EXT", lgth))
)

##### 2. Plot analysis #####

apmPlot <- ggplot(apm, aes(x=Sequence, fill = Member)) +
  geom_bar(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0, lgth) +
  labs(x = "States",
       title = "APM")

amPlot <- ggplot(am, aes(x=Sequence, fill = Member)) +
  geom_bar(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0, lgth) +
  labs(x = "States",
       title = "AM")

pmPlot <- ggplot(pm, aes(x=Sequence, fill = Member)) +
  geom_bar(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0, lgth) +
  labs(x = "States",
       title = "PM")

equPlot <- ggplot(equ, aes(x=Sequence, fill = Member)) +
  geom_bar(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0, lgth) +
  labs(x = "States",
       title = "EQU")

extPlot <- ggplot(ext, aes(x=Sequence, fill = Member)) +
  geom_bar(show.legend = FALSE) +
  facet_wrap(~Member, scales = "free") +
  ylim(0, lgth) +
  labs(x = "States",
       title = "EXT")

mainBox <- textGrob("States histogram",gp=gpar(fontsize=20,font=3))
grid.arrange(apmPlot, amPlot, pmPlot, equPlot, extPlot, ncol=3, nrow=2, top = mainBox)
#dev.print(pdf, file = "StatesHist.pdf") 

dev.print(png, file = "StatesHist.png", width = 1000, height = 800) 


##### 2. GOF #####

members <- c("FM", "SM", "Global")
memNames <- members 

distance <- matrix(NA, nrow = length(models), ncol = length(members), byrow = TRUE)
pvalue <- matrix(NA, nrow = length(models), ncol = length(members), byrow = TRUE)


for(j in 1:length(models)){
  empX <- countEmp(states = s, chainFM = sequence[,1,j], chainSM = sequence[,2,j])
  empY <- countEmp(states = s, chainFM = sequence[,2,j], chainSM = sequence[,1,j])
  theoX <- countTheo(lambda = lambda[,,j,1], countEmp = empX, condition = "U")
  theoY <- countTheo(lambda = lambda[,,j,2], countEmp = empY, condition = "U")
  compX <- accuracyLocal(population = theoX, empirical = empX)
  compY <- accuracyLocal(population = theoY, empirical = empY)
  comp <- accuracyGlobal(populationFM = theoX, empiricalFM = empX, populationSM = theoY, empiricalSM = empY)
  dist <- c(compX$statistic[["X-squared"]], compY$statistic[["X-squared"]], comp$statistic[["X-squared"]])
  value <- c(compX$p.value, compY$p.value, comp$p.value)
  distance[j,] <- dist
  pvalue[j,] <- value 
}

assign("distanceGOF", distance)
assign("pvalueGOF", pvalue)

lambdaX <- lambda[,,2,1]
lambdaY <- lambda[,,2,2]
countEmpX <- countEmp(states = 2, chainFM = as.numeric(am[1:1000,1]), chainSM = as.numeric(am[1001:2000,1]))

countTheoGOF <- countTheo(lambda = lambdaX, countEmp = countEmpX, condition = "U", type = "AM")
countTheoLRT <- countTheo(lambda = lambdaX, countEmp = countEmpX, condition = "C", type = "AM")

round(countTheoGOF, 3)
round(countTheoLRT, 3)

############
### Save ###
############

save(two100, file = "two100.Rdata")
save(two500, file = "two500.Rdata")
save(two1000, file = "two1000.Rdata")
save(two10000, file = "two10000.Rdata")

save(three100, file = "three100.Rdata")
save(three500, file = "three500.Rdata")
save(three1000, file = "three1000.Rdata")
save(three10000, file = "three10000.Rdata")

save(twodyadAPM, file = "twodyadAPM.Rdata")
save(twodyadPATTERN, file = "twodyadPATTERN.Rdata")
save(twodyadEQU, file = "twodyadEQU.Rdata")
save(twodyadEXT, file = "twodyadEXT.Rdata")

save(threedyadAPM, file = "threedyadAPM.Rdata")
save(threedyadPATTERN, file = "threedyadPATTERN.Rdata")
save(threedyadEQU, file = "threedyadEQU.Rdata")
save(threedyadEXT, file = "threedyadEXT.Rdata")

save(tableGOFAPM, file = "tableGOFAPM.Rdata")
save(tableGOFEQU, file = "tableGOFEQU.Rdata")
save(tableGOFEXT, file = "tableGOFEXT.Rdata")

save(testGOFAPM, file = "testGOFAPM.Rdata")
save(testGOFAPM, file = "testGOFAPM.Rdata")
save(testGOFAPM, file = "testGOFAPM.Rdata")

save(lambda, file = "empiricalLambda.Rdata")
save(apm, file = "empiricalAPM.Rdata")
save(am, file = "empiricalAM.Rdata")
save(pm, file = "empiricalPM.Rdata")
save(equ, file = "empiricalEQU.Rdata")
save(ext, file = "empiricalEXT.Rdata")


save(distanceGOF, file = "distanceGOF.Rdata")
save(pvalueGOF, file = "pvalueGOF.Rdata")

