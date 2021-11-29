############ THR model #########################

##### (1) Set-up Model ######
### PACKAGES USED
if(!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

if(!require(reshape2)) install.packages('reshape2')
library(reshape2)

set.seed(1234)

# Reading the data needed from csv files
life.tables <- read.csv("life-table.csv", header=TRUE)         # age group, sex
hazard.function <- read.csv("hazardfunction.csv", header=TRUE) # from the regression analysis: lngamma,cons,age,male,NP1
cov_coeff <- read.csv("cov55.csv", row.names=1, header=TRUE)      # explanatory variable covariance matrix

### Structural inputs
state.names <- c("PrimaryTHR","SuccessP","RevisionTHR","SuccessR","Death")
n.states <- length(state.names)         # number of states in the model
seed <- c(1,0,0,0,0)                    # starting states of the model
cycles <- 60                            # number of cycles running the model
cycle.v <- 1:cycles                     # a vector of cycle numbers 1 - 60
cDR <- 0.06                             # set the discount rate for costs (6%)
oDR <- 0.015                            # set the discount rate for outcomes (15%)
discount.factor.c <- 1/(1+cDR)^cycle.v  # the discount factor matrix
discount.factor.o <- 1/(1+oDR)^cycle.v  # discount factor matrix for utility 

##### (2) DETERMINISTIC PARAMETERS ######

# COSTS
cSP0 <- 394    # standard prosthesis
cNP1 <- 579    # new prosthesis 1
cPrimary <- 0  # primary THR procedure - set to 0 for this model 
cSuccess <- 0  # success - set to 0 for this model


### (3) Define Hazard function ####

## Coefficients - on the log hazard scale
# vector of mean values from the regression analysis
haz_coeff <- hazard.function$coefficient


#### (4) Defining shape and scale parameters ####
##NG: I'm not sure this is what they're called in a Beta distn

beta_params <- list()

# revision surgery
mn.cRevision <- 5294 # mean cost
se.cRevision <- 1487 # standard error of cost

# translate for ~Beta(a, b)
a.cRevision <- (mn.cRevision/se.cRevision)^2  # alpha
b.cRevision <- (se.cRevision^2)/mn.cRevision  # beta

beta_params$cRevision <- c(a = a.cRevision, b = b.cRevision)

### Transition probabilities

# operative mortality from primary surgery
a.PTHR2dead <- 2                  # alpha
b.PTHR2dead <- 100 - a.PTHR2dead  # beta

beta_params$PTHR2dead <- c(a = a.PTHR2dead, b = b.PTHR2dead)

# re-revision risk
a.rrr <- 4            # alpha 
b.rrr <- 100 - a.rrr  # beta

beta_params$rrr <- c(a = a.rrr, b = b.rrr)

## UTILITIES
# primary prosthesis
mn.uSuccessP <- 0.85                                             # mean utility value
se.uSuccessP <- 0.03                                             # standard error utility value
ab.uSuccessP <- mn.uSuccessP*(1-mn.uSuccessP)/(se.uSuccessP^2)-1 # estimating alpha plus beta (ab)
a.uSuccessP <- mn.uSuccessP*ab.uSuccessP                         # estimating alpha (a)
b.uSuccessP <- a.uSuccessP*(1-mn.uSuccessP)/mn.uSuccessP         # estimating beta (b)

beta_params$uSuccessP <- c(a = a.uSuccessP, b = b.uSuccessP, ab = ab.uSuccessP)

## revision surgery
mn.uSuccessR <- 0.75                                             # mean utility value
se.uSuccessR <- 0.04                                             # standard error utility value
ab.uSuccessR <- mn.uSuccessR*(1-mn.uSuccessR)/(se.uSuccessR^2)-1 # alpha + beta (ab)
a.uSuccessR <- mn.uSuccessR*ab.uSuccessR                         # alpha (a)
b.uSuccessR <- a.uSuccessR*(1-mn.uSuccessR)/mn.uSuccessR         # beta(b)

beta_params$uSuccessR <- c(a = a.uSuccessR, b = b.uSuccessR, ab = ab.uSuccessR)

## during the revision period 
mn.uRevision <- 0.30                                             # mean utility score
se.uRevision <- 0.03                                             # standard error utility score
ab.uRevision <- mn.uRevision*(1-mn.uRevision)/(se.uRevision^2)-1 # alpha + beta (ab)
a.uRevision  <- mn.uRevision*ab.uRevision                        # alpha (a)
b.uRevision  <- a.uRevision*(1-mn.uRevision)/mn.uRevision        # beta(b)

beta_params$uRevision <- c(a = a.uRevision, b = b.uRevision, ab = ab.uRevision)

beta_params



#####**** (5) SET UP SAMPLE FUNCTION ******#####

sim.runs <- 100

psa.sampling <- function(age = 60,
                         male = 0, #NG: why not FALSE?
                         beta_params,
                         haz_coeff,
                         cov_coeff,
                         sim.runs = 1000) {
  
  #### description: sample probabilistic parameters according to age and sex
  #### params: age (numeric), male (0 for female, 1 for male), though 
  ####         other variables which are defined above are called within the function
  ####         e.g. cycles
  #### return: a list with data frames and vectors for probabilistic parameters 
  
  ###  Transition probabilities
  ##NG: not a good idea to implicitly refer to enclosing environment
  tp.PTHR2dead <- rbeta(sim.runs, beta_params$PTHR2dead['a'], beta_params$PTHR2dead['b']) ## OMR following primary THR
  
  ##NG: should this be the same as above?
  tp.RTHR2dead <- rbeta(sim.runs, beta_params$PTHR2dead['a'], beta_params$PTHR2dead['b'])  ## OMR following revision THR
  
  ## creating a data.frame of sampled transition probabilities
  omr.df <- data.frame(tp.PTHR2dead, tp.RTHR2dead)
  
  tp.rrr.vec <- rbeta(sim.runs, beta_params$rrr['a'], beta_params$rrr['b']) ## Re-revision risk transitions vector
  
  ###  Costs
  c.revision.vec <- rgamma(sim.runs,
                           shape = beta_params$cRevision['a'],
                           scale = beta_params$cRevision['b']) ## Gamma distribution draw for cost of revision surgery
  
  ## Utilities
  uSuccessP <- rbeta(sim.runs, beta_params$uSuccessP['a'], beta_params$uSuccessP['b']) 
  uSuccessR <- rbeta(sim.runs, beta_params$uSuccessR['a'], beta_params$uSuccessR['b']) 
  uRevision <- rbeta(sim.runs, beta_params$uRevision['a'], beta_params$uRevision['b']) 
  
  ## Make a data frame to pass into the function
  state.utilities.df <- data.frame(uprimary = rep(0, sim.runs), 
                                   uSuccessP,
                                   uRevision,
                                   uSuccessR,
                                   udeath = rep(0, sim.runs))
  
  ## (6) Hazard function sampling ####
  n_draws <- 5
  z <-
    matrix(
      rnorm(n_draws*sim.runs, 0, 1),                   # normal distribution N(0,1)
      nrow = sim.runs, ncol = n_draws)                 # n_draws random samples, by sim.runs
  r.table <- matrix(0, nrow = sim.runs, ncol = n_draws)
  colnames(r.table) <- names(cov_coeff)
  
  cholm <- t(chol(t(cov_coeff)))             # lower triangle of the Cholesky decomposition
  
  for (i in 1:sim.runs) {
    Tz <- cholm %*% z[i, ] 
    x <- haz_coeff + Tz 
    r.table[i, ] <- x[, 1]
  }
  
  r <- as.data.frame(r.table)
  gamma.vec <- exp(r$lngamma)
  lambda.vec <- exp(r$cons + age * r$age + male*r$male)
  RR.vec <- exp(r$NP1)
  survival.df <- data.frame(gamma.vec, lambda.vec)
  
  #### (7) Life-table sampling #####
  
  # set life table values based on age and sex (not probabilistic but dependent on)
  # age and sex variables chosen
  current.age <- age + cycle.v ## a vector of cohort age throughout the model
  interval <- findInterval(current.age, life.tables$Index)
  death.risk <- data.frame(age = current.age,
                           Males = life.tables[interval, 'Males'],
                           Females = life.tables[interval, 'Females'])
  col.key <- ifelse(male, "Males", "Females")
  mortality.vec <- death.risk[, col.key]
  
  ## combine outputs
  list(RR.vec = RR.vec,
       omr.df = omr.df,
       tp.rrr.vec = tp.rrr.vec,
       survival.df = survival.df,
       c.revision.vec = c.revision.vec,
       state.utilities.df = state.utilities.df,
       mortality.vec = mortality.vec)
}


sample.output <-
  psa.sampling(beta_params = beta_params,
               haz_coeff = haz_coeff,
               cov_coeff = cov_coeff,
               sim.runs = sim.runs)


##NG: this is as far as I've got (12/11/2021)#





### (8) Defining outputs from sampling #####
RR.vec <- sample.output$RR.vec
omr.df <- sample.output$omr.df
tp.rrr.vec <- sample.output$tp.rrr.vec
survival.df <- sample.output$survival.df
c.revision.vec <- sample.output$c.revision.vec 
state.utilities.df <- sample.output$state.utilities.df 
mortality.vec <- sample.output$mortality.vec
## take a look at the above vectors and data.frames to see the samples produced
## these are defined out of the list here for use in the value of information analysis later on

####***** (9) THR MODEL FUNCTION ****#####
model.THR <- function(RR.NP1,      # from RR.vec
                      omr,         # from omr.df
                      tp.rrr,      # from tp.rrr.vec
                      survival,    # from survival.df
                      c.revision,  # from c.revision.vec
                      state.util,  # from state.utilities.df
                      mortality.vec) { ## background mortality based on life tables for age/sex combination
  ### FUNCTION: Run the THR model on the deterministic and probablistic parameters
  ## INPUTS: 3 data.frame rows (1 row from omr.df, survival.df & state.utilities.df) 
  ##        and 3 vector values (1 value from vectors; c.revision.vector, tp.rrr.vector, RR.vec)though 
  #### other variables which are defined above are called within the function
  ### e.g. cycles & age/sex are set in the psa.sampling function
  ## OUTPUTS: a labelled numeric output providing costs and qalys for SP0 and NP1, given the set of input parameters
  
  ## First, we need to unpack values from the data.frames provided to the function 
  # COSTS:
  state.costs<-c(cPrimary, cSuccess, c.revision, cSuccess, 0) ## a vector with the costs for each state
  # UTILITIES:
  state.utilities <- unlist(state.util)
  # TRANSITIONS
  tp.PTHR2dead <- unlist(omr[1,1])
  tp.RTHR2dead <- unlist(omr[1,2])
  gamma <- unlist(survival[1,1])
  lambda <- unlist(survival[1,2])
  
  ### (10) Integrating revision and mortality risks #### 
  revision.risk.sp0 <- 1- exp(lambda * ((cycle.v-1) ^gamma-cycle.v ^gamma))
  revision.risk.np1 <- 1- exp(lambda * RR.NP1 * ((cycle.v-1) ^gamma-cycle.v ^gamma))
  # Transition arrays
  tm.SP0 <- array(data=0,dim=c(n.states, n.states, cycles),
                  dimnames= list(state.names, state.names, 1:cycles)) ## an empty array of dimenions (number of states, number of states, number of cycles)
  # Now using vectorisation to complete this 
  tm.SP0["PrimaryTHR","Death",] <- tp.PTHR2dead 
  tm.SP0["PrimaryTHR","SuccessP",] <- 1 - tp.PTHR2dead 
  tm.SP0["SuccessP","RevisionTHR",] <- revision.risk.sp0 
  tm.SP0["SuccessP","Death",] <- mortality.vec
  tm.SP0["SuccessP","SuccessP",] <- 1 - revision.risk.sp0 - mortality.vec
  tm.SP0["RevisionTHR","Death",] <- tp.RTHR2dead + mortality.vec
  tm.SP0["RevisionTHR","SuccessR",] <- 1 - tp.RTHR2dead - mortality.vec
  tm.SP0["SuccessR","RevisionTHR",] <- tp.rrr
  tm.SP0["SuccessR","Death",] <- mortality.vec
  tm.SP0["SuccessR","SuccessR",] <- 1 - tp.rrr - mortality.vec
  tm.SP0["Death","Death",] <- 1 
  #  Create a trace for the standard prosthesis arm
  trace.SP0 <- matrix(data=0, nrow=cycles, ncol=n.states)
  colnames(trace.SP0) <- state.names
  
  trace.SP0[1,] <- seed%*%tm.SP0[,,1]
  
  for (i in 2:cycles) trace.SP0[i,] <- trace.SP0[i-1,]%*%tm.SP0[,,i]
  #### (11)  NP1 ARM #####
  tm.NP1 <- array(data=0,dim=c(n.states, n.states, cycles),
                  dimnames= list(state.names, state.names, 1:cycles)) ## an empty array of dimenions (number of states, number of states, number of cycles)
  ## naming all dimensions
  ### create a loop that creates a time dependent transition matrix for each cycle
  tm.NP1["PrimaryTHR","Death",] <- tp.PTHR2dead ## Primary THR either enter the death state or.. or..
  tm.NP1["PrimaryTHR","SuccessP",] <- 1 - tp.PTHR2dead ## they go into the success THR state 
  tm.NP1["SuccessP","RevisionTHR",] <- revision.risk.np1 ## revision risk with NP1 treatment arm 
  tm.NP1["SuccessP","Death",] <- mortality.vec
  tm.NP1["SuccessP","SuccessP",] <- 1 - revision.risk.np1 - mortality.vec
  tm.NP1["RevisionTHR","Death",] <- tp.RTHR2dead + mortality.vec
  tm.NP1["RevisionTHR","SuccessR",] <- 1 - tp.RTHR2dead - mortality.vec
  tm.NP1["SuccessR","RevisionTHR",] <- tp.rrr
  tm.NP1["SuccessR","Death",] <- mortality.vec[i]
  tm.NP1["SuccessR","SuccessR",] <- 1 - tp.rrr - mortality.vec
  tm.NP1["Death","Death",] <- 1 
  
  #  Create a trace for the standard prosthesis arm
  trace.NP1 <- matrix(data=0,nrow=cycles,ncol=n.states)
  colnames(trace.NP1) <- state.names
  trace.NP1[1,] <- seed%*%tm.NP1[,,1]
  for (i in 2:cycles) trace.NP1[i,] <- trace.NP1[i-1,]%*%tm.NP1[,,i]
  
  #### (12) Analysis #####
  # COST #
  cost.SP0 <- trace.SP0%*%state.costs  
  disc.cost.SP0 <- (discount.factor.c%*%cost.SP0) + cSP0   
  cost.NP1 <- trace.NP1%*%state.costs  
  disc.cost.NP1 <- (discount.factor.c%*%cost.NP1) + cNP1 
  
  # QALYS #
  QALYs.SP0 <- trace.SP0%*%state.utilities ## utility per cycle
  disc.QALYs.SP0 <- colSums(discount.factor.o%*%QALYs.SP0) ## total discounted utility
  QALYs.NP1 <- trace.NP1%*%state.utilities ## utility per cycle
  disc.QALYs.NP1 <- colSums(discount.factor.o%*%QALYs.NP1) ## total discounted utility
  
  c(cost.SP0 = disc.cost.SP0,
    qalys.SP0 = disc.QALYs.SP0,
    cost.NP1 = disc.cost.NP1,
    qalys.NP1 = disc.QALYs.NP1,
    inc.cost = disc.cost.NP1 - disc.cost.SP0,
    inc.qalys = disc.QALYs.NP1 - disc.QALYs.SP0)
}

#### (13) Running the simulations ########

## creating an empty data.frame for simulation results to fill:
simulation.results <- data.frame("cost.SP0" = rep(as.numeric(NA), sim.runs), ## rep() creates sim.runs rows of values
                                 "qalys.SP0"= rep(as.numeric(NA), sim.runs),
                                 "cost.NP1" = rep(as.numeric(NA), sim.runs),
                                 "qalys.NP1" = rep(as.numeric(NA), sim.runs),
                                 "inc.cost" = rep(as.numeric(NA), sim.runs),
                                 "inc.qalys"=  rep(as.numeric(NA), sim.runs))


## running the simulations and filling the simulation.results data.frame:
for(i in 1:sim.runs){
  simulation.results[i,] <- model.THR(RR.vec[i],
                                      omr.df[i,],  
                                      tp.rrr.vec[i],
                                      survival.df[i,],
                                      c.revision.vec[i], 
                                      state.utilities.df[i,],
                                      mortality.vec) 
}


### (14) Estimating average net monetary benefit ####

p.CE <- function(WTP, simulation.results) {# Estimate the probability of cost-effectiveness for a given willingness-to-pay ceiling ratio
  ## a function that estimates the probability of the new intervention
  # being cost-effective 
  # INPUTS: WTP = willingness to pay value (numeric)
  #         simulation.results = a data.frame output of PSA simulations which includes
  #         columns names "inc.qalys" and "inc.cost"   
  # OUTPUTS: A numeric value specifying the probability of cost-effectiveness given the inputs
  nmb <- simulation.results[,"inc.qalys"]*WTP - simulation.results[,"inc.cost"] ## vector of NMB estimates for each simulation length 1:sim.runs
  CE <- nmb > 0   ## vector of TRUE/FALSE for each simulation length 1:sim.runs
  probCE <- mean(CE) ## the mean value of TRUE (=1) and FALSE (=0)
  return(probCE)
}

# Generate CEAC table
WTP.values <- seq(from = 0, to = 50000, by = 10) ## use the seq() function to get a vector of specified numeric values
CEAC <- data.frame(WTP = WTP.values, 
                   pCE = rep(as.numeric(NA),length(WTP.values)))

for (i in 1:length(WTP.values)) {
  CEAC[i,"pCE"] <- p.CE(CEAC[i,"WTP"], simulation.results)
}

##### (15) Subgroup Analyses ######

# CREATE ARRAY TO STORE THE RESULTS OF THE MODEL IN EACH SUBGROUP
subgroups.names <- c("Male 40", "Male 60", "Male 80", "Female 40", "Female 60", "Female 80")
subgroups.n <- length(subgroups.names)
simulation.subgroups <- array(data = 0,
                              dim = c(sim.runs, length(colnames(simulation.results)), subgroups.n),
                              dimnames = list(1:sim.runs, colnames(simulation.results), subgroups.names))

# Run model for each subgroup, inputting the age and sex into the function, and record results within the array
sample.sub <- list()
sample.sub[[1]] <- psa.sampling(age = 40, male = 1)
sample.sub[[2]] <- psa.sampling(age = 60, male = 1)
sample.sub[[3]] <- psa.sampling(age = 80, male = 1)
sample.sub[[4]] <- psa.sampling(age = 40, male = 0)
sample.sub[[5]] <- psa.sampling(age = 60, male = 0)
sample.sub[[6]] <- psa.sampling(age = 80, male = 0)

for(i in 1:sim.runs){ 
  for(j in 1:6){ ## column = each subgroup
    simulation.subgroups[i,,j] <- model.THR(sample.sub[[j]]$RR.vec[i],
                                            sample.sub[[j]]$omr.df[i,],
                                            sample.sub[[j]]$tp.rrr.vec[i], 
                                            sample.sub[[j]]$survival.df[i,],
                                            sample.sub[[j]]$c.revision.vec[i], 
                                            sample.sub[[j]]$state.utilities.df[i,],
                                            mortality.vec = sample.sub[[j]]$mortality.vec)
  }
}

# Create a CEAC table with lambda value sequence
WTP.values <- seq(from = 0, to = 50000, by = 50)
CEAC.subgroups <- matrix(data= as.numeric(NA), nrow=length(WTP.values), ncol=subgroups.n + 1)
CEAC.subgroups <- as.data.frame(CEAC.subgroups)
colnames(CEAC.subgroups) <- c("WTP", subgroups.names)
# Estimate probability cost-effective for all subgroups
for (i in 1:length(WTP.values)) {
  CEAC.subgroups[i,1] <- WTP.values[i]
  CEAC.subgroups[i,2] <- p.CE(WTP.values[i], simulation.subgroups[,,1])
  CEAC.subgroups[i,3] <- p.CE(WTP.values[i], simulation.subgroups[,,2])
  CEAC.subgroups[i,4] <- p.CE(WTP.values[i], simulation.subgroups[,,3])
  CEAC.subgroups[i,5] <- p.CE(WTP.values[i], simulation.subgroups[,,4])
  CEAC.subgroups[i,6] <- p.CE(WTP.values[i], simulation.subgroups[,,5])
  CEAC.subgroups[i,7] <- p.CE(WTP.values[i], simulation.subgroups[,,6])
}


######***PLOTS****#####################

#### (16) COST-EFFECTIVENESS PLANE #####
## Plotting:
xlabel = "Incremental QALYs"
ylabel = "Incremental costs"
ggplot(simulation.results) + 
  geom_point(shape = 21, size = 2, colour = "black", fill = NA, alpha = 0.5, aes(x=inc.qalys, y=inc.cost)) + 
  labs(x = xlabel, text = element_text(size=10)) + labs (y = ylabel, text = element_text(size=10)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12),
        plot.margin=unit(c(1.2,0.5,0,1.2),"cm"))


## plotting:
xlabel = "Willingness to pay threshold"
ylabel = "Probability cost-effective"
ggplot(CEAC) + geom_line(aes(x=WTP, y=pCE), size=1) + 
  labs(x = xlabel, text = element_text(size=10)) + labs(y = ylabel, text = element_text(size=10)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12)) + 
  scale_x_continuous(expand = c(0, 0.1)) + 
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0))

### (17) COST-EFFECTIVENESS ACCEPTABILITY CURVES ####
## To plot the CEAC for subgroups We need to reshape the data from wide to long to use in ggplot 
CEAC.subgroups.long <- melt(CEAC.subgroups, id.vars = c("WTP"))
colnames(CEAC.subgroups.long) <- c("WTP", "group", "pCE")
## plotting:
xlabel = "Willingness to pay threshold"
ylabel = "Probability cost-effective"
ggplot(CEAC.subgroups.long) + geom_line(aes(x=WTP, y=pCE, color=group), size=1) + 
  labs(x = xlabel, text = element_text(size=10)) + labs(y = ylabel, text = element_text(size=10)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12)) + 
  scale_x_continuous(expand = c(0, 0.1)) + 
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0))
