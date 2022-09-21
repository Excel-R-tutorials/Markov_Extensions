########### THR MODEL: VALUE OF INFORMATION ##############

##### (1) Set up model script and source THR model functions  #######

# Loading in data and model
# this also loads libraries called through "THR_Model.R"
source("THR_Model.R")
## Note that simulation results are read in from this, and used throughout 
# see head(simulation.results)

#### (2) Setting value of information population parameters ####

population <- 40000 
years <- 10
evpi.disc <- 0.06
population.seq <- population * (1/(1+evpi.disc) ^ c(0:(years-1)))
effective.population <- sum(population.seq)
## Create a vector of willingness to pay values
WTP.values <- seq(from = 0, to = 50000, by = 50)


#### (3) Expected value of perfect information (EVPI) at the population level ####

est.EVPI.pop <-function(WTP, effective.population, simulation.results) {

  # Estimate the NMB for two treatments, for each simulation
  nmb.SP0 <- simulation.results$qalys.SP0 * WTP - simulation.results$cost.SP0
  nmb.NP1 <- simulation.results$qalys.NP1 * WTP - simulation.results$cost.NP1
  nmb.table <- data.frame(nmb.SP0, nmb.NP1)
  
  # Estimating NMB with current information
  av.current <- apply(nmb.table, 2, mean)
  max.current <- max(av.current)
  
  # Estimating NMB with perfect information
  max.nmb <- apply(nmb.table, 1, max) 
  av.perfect <- mean(max.nmb)
  
  # Generating EVPI values
  EVPI.indiv <- av.perfect - max.current
  pop.EVPI <- effective.population * EVPI.indiv
  
  return(pop.EVPI)
  
}

## create a data frame to capture the WTP values and subsequent population EVPI results 

EVPI.results <- data.frame(WTP=WTP.values, EVPI=NA)

for (i in 1:length(WTP.values)) {
  EVPI.results[i,2] <- est.EVPI.pop(WTP.values[i], effective.population, simulation.results)
}

## Show the highest EVPI value
EVPI.results[EVPI.results$EVPI == max(EVPI.results$EVPI),] 

#### (4) Setting up the EVPPI inner and outer loop framework #### 

## Enter inner and outer loop numbers 
inner.loops <- 100
outer.loops <- 100

# Create a matrix to store the inner loop results
inner.results <- matrix(0, inner.loops, 6)
colnames(inner.results) <- c("Cost SP0", "QALY SP0", "Cost NP1", "QALY NP1", "Inc Cost", "Inc QALY")

evppi.results.SP0 <- matrix(0, nrow = outer.loops, ncol = length(WTP.values))
colnames(evppi.results.SP0) <- as.character(WTP.values)
evppi.results.NP1 <- evppi.results.SP0 

##### (5) The Net Monetary Benefit Function #######

# Creating a function to estimate NMB across WTP values (for inner loop results)
nmb.function <- function(WTP, results){

  nmb.table <- matrix(WTP, ncol = length(WTP), nrow = dim(results)[1],  byrow = TRUE) 
  
  SP0 <- ((results[,"QALY SP0"] * nmb.table) - results[,"Cost SP0"])  
  NP1 <- ((results[,"QALY NP1"] * nmb.table) - results[,"Cost NP1"])
  
  nmb.SP0 <- apply(SP0, 2, mean)
  nmb.NP1 <- apply(NP1, 2, mean) 
    
  ###  OUTPUTS: a list of the NMB under SP0 (nmb.SP0) and NP1 (nmb.NP1)
  return(list(nmb.SP0, nmb.NP1))
  
}

##### (6) Calculating EVPPI values for Willingness to Pay Values  #####

## Function to estimate EVPPI across WTP values, using the outer loop results
gen.evppi.results <- function(evppi.results1 = evppi.results.SP0, evppi.results2 = evppi.results.NP1, WTP = WTP.values, effective.pop = effective.population){

  ## calculate the mean NMB for current and new treatments, at each WTP 
  current.info1 <- apply(evppi.results1, 2, mean)
  current.info2 <- apply(evppi.results2, 2, mean)
  
  ## calculate the max NMB (overall) from either current or new treatments, at each WTP
  current.info <- pmax(current.info1, current.info2)
  
  # Create an array 
  evppi.array <- array(0, dim = c(outer.loops, length(WTP), 2))
  evppi.array[,,1] <- evppi.results1
  evppi.array[,,2] <- evppi.results2
  
  # calculate the max NMB for each treatment, at each WTP, for each simulation 
  # This is so that the best treatment can be select for each simulation
  perf.info.sims <- apply(evppi.array, c(1,2), max)
  perf.info <- apply(perf.info.sims, 2, mean)
  
  # Work out the difference between perfect information (for each simulation) vs. current information, at each WTP
  evppi.results <- c(perf.info - current.info) * effective.pop
  
  ### OUTPUT: EVPPI results     
   return(evppi.results)
  
}

#### (7) Run all EVPPI simulations  ######

## Now the EVPPI loops will be run - each selected different values for inner and outer loops
parameter.groups <- 6 ## number of parameters/parameter groups to run the EVPPI on
evppi.wide <- data.frame(wtp = WTP.values, RR = NA, OMR = NA,  surv = NA, r.cost = NA, RRR = NA, util = NA)
colnames(evppi.wide) <- c("WTP", "NP1 Relative risk",  "Operative mortality ratios", "Survival parameters", "Revision cost", "Re-revision risk", "Utilities")

# Set progres bar 
pb = txtProgressBar(min = 0, max = outer.loops*parameter.groups, initial = 0, style = 3)

# Create loops for 1) parameter groups, 2) outer loops 3) inner loops
for(j in 1:parameter.groups){
  
  for(a in 1:outer.loops){
    
    for(b in 1:inner.loops){
      
      # The 'partial' parameter will be included in the outer loop ('a')
      # The other parameters are included in the inner loop ('b')
      
      if(j==1) rr.n <- a else rr.n <- b
      if(j==2) omr.n <- a else omr.n <- b 
      if(j==3) surv.n <- a else surv.n <- b   
      if(j==4) c.rev.n <- a else c.rev.n <- b      
      if(j==5) rrr.n <- a else rrr.n <- b
      if(j==6) util.n <- a else util.n <- b
      
      inner.results[b,] <-  model.THR(RR.vec[rr.n], omr.df[omr.n,],  tp.rrr.vec[rrr.n], 
                                      survival.df[surv.n,],c.revision.vec[c.rev.n], 
                                      state.utilities.df[util.n,], mortality.vec = mortality.vec) 
      
    }
    
    #after each inner loop PSA, calculate the mean NMB for each treatment compactor and store the results
    nmb <- nmb.function(WTP.values, inner.results)
    
    evppi.results.SP0[a,] <- nmb[[1]]
    evppi.results.NP1[a,] <- nmb[[2]]
    setTxtProgressBar(pb,(j-1)*outer.loops+a)
    
  }
  
  evppi.wide[,j+1] <- gen.evppi.results() 
  
}

# Preview of results 
head(evppi.wide,20)

######*** Graphical Outputs ****#####################

##### (a) Ploting results EVPI, per population ####
ggplot(EVPI.results) + geom_line(aes(x=WTP, y=EVPI), size=1) + 
  labs(x = "Willingness to pay threshold", text = element_text(size=10)) + 
  labs(y = "Expected Value of Perfect Information", text = element_text(size=10)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
        plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
  scale_x_continuous(labels = scales::comma, expand = c(0, 0.1)) + 
  scale_y_continuous(labels = scales::comma, expand = c(0, 0))

#### (b) PLOTTING EVPPI RESULTS ACROSS WILLINGNESS-TO-PAY THRESHOLDS ####

# Convert from wide to long format
evppi.long <- reshape2::melt(evppi.wide, id.vars = c("WTP"))

ggplot(evppi.long) + geom_line(aes(x=WTP, y=value, colour = variable), size=0.75) + 
  labs(x = "Willingness to pay threshold", text = element_text(size=10)) + 
  labs(y = "EVPPI", text = element_text(size=10)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12),
        plot.margin=unit(c(0.5,0,0,0.5),"cm")) + 
  scale_x_continuous(labels = scales::comma, limits = c(0, 10000), expand = c(0, 0.1)) + 
  scale_y_continuous(labels = scales::comma, expand = c(0, 0))
# Note you will get a warning here that the plot does not show all the data in the data.frame



#### (c) PLOTTING EVPPI RESULTS AT A WILLINGNESS-TO-PAY THRESHOLD VALUE ####

sub.evppi <- subset(evppi.long, WTP==2200)

ggplot(sub.evppi, aes(x=variable, y=value)) +
  geom_bar(stat="identity", fill="steelblue")+
  labs(x = "Parameter Group", text = element_text(size=4)) + 
  labs(y = "EVPPI", text = element_text(size=4)) + theme_classic() +
  theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)),
        axis.text.x=element_text(angle=45,hjust=1), 
        panel.grid.major = element_line(), panel.grid.minor = element_line(), 
        legend.key.width=unit(1.8,"line"), text = element_text(size=12)) + 
  scale_y_continuous(labels = scales::comma, expand = c(0, 0))
