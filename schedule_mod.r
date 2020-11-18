library(tidyverse)
library(grid)

##-------Problem and Algorithm Parameter Definition Section------------------

                           
# Number of supervisors available
nsupervisorsP <- seq(1,10,4)

# Number of non-supervisor employees
nemployeesP <- seq(1,30, 8)

# Number of days in the schedule
ndaysP <- c(7, 14, 30)

# TRUE/FALSE: does the schedule repeat? (so for example working
# last night into first morning is penalized)
periodic <- TRUE

# Number of staff present for each shift
# e.g., c(2, 3, 3) has three shifts with 2, 3, 3 people respectively
shiftstaffingP <- list(c(2,3), c(2,2,2,2,2), c(10,10,10))

# Penalties for violating each of our four constraints
# Penalty 1: Each shift without a supervisor present
# Penalty 2: Each instance of consecutive shifts for same person (including overnight)
# Penalty 3: Each instance of three days worked in a row
# Penalty 4: Weight for the stdev of each employee's shifts in the schedule
#            (trying to make everyone li similar amounts)
penaltiesP <- list(c(1, 1, 1, 1), c(1, 1, 1, 3), c(1, 5, 5, 1))

# Number of iterations (if too small, might not be able to explore enough)
niterP <- c(5e+4, 1e+5, 5e+5)

# Different types of cooling schedule, details listed in the next section
coolingP <- c("greedy", "lin", "log", "expo")

# Different step sizes - the # of staff on the schedules to be changed
stepP <- seq(1,10,3)


# Generate unique combinations of the parameters listed above  
parameter_df <- expand.grid(nsupervisors = nsupervisorsP, 
                            nemployees = nemployeesP, 
                            ndays = ndaysP, 
                            shiftstaffing = shiftstaffingP,
                            penalties = penaltiesP,
                            niter = niterP,
                            cooling = coolingP,
                            step = stepP)

# If total amount of shifts per day is more than total staff, then it is not possible
sc <- rep(0, dim(parameter_df)[1])
for (i in 1:dim(parameter_df)[1]) {
  sc[i] <- sum(Reduce(`+`, parameter_df$shiftstaffing[i]))
}

parameter_df$shiftsSum <- sc

parameter_df <- parameter_df %>% 
  filter(shiftsSum < (nsupervisors + nemployees)) %>% 
  select(-contains("shiftsSum"))

# Change data type
parameter_df$cooling <- as.character(parameter_df$cooling)

# Prepare the dataframe to store more incoming data 
parameter_df$energy <- NA
parameter_df$min_energy <- NA
parameter_df$first_min_R <- NA
parameter_df$final_energy <- NA
parameter_df$acceptsR <- NA

parameter_df$Index <- as.numeric(rownames(parameter_df))
# Decrease the total instances (high runtime)
# 
# parameter_df <- parameter_df[sample(nrow(parameter_df), (dim(parameter_df)[1] / 10)), ]
# row.names(parameter_df) <- NULL
## END OF PARAMETER SECTION


##-------Function Declaration Section---------------------------

# Function to run the algorithm based on the prepared parameter dataframe

simuated_annealing <- function(index)
{ 
  # The cooling schedule: what should the temperature be at step i out of n
  # 
  coolingschedule <- function(i, n, method)
  {
    result = switch(
      method,
      "greedy" = 1e-3,
      "lin" = 1 - (i/(n+1)),
      "log" = 0.02/log(1+(i+1)/n),
      "expo" = n^((1-i)/n)
    )
    return(result)
  }
  
  # Function to calculate the "energy" (total cost of violations)
  get_energy <- function(x, penalties)
  { 
    energy <- 0
    for (i in 1:ndays)
    {
      for (j in 1:length(shiftstaffing))
      {
        # No supervisor present on the shift
        if (sum(x[shiftindices[j]:(shiftindices[j+1]-1),i]<=nsupervisors) == 0)
        {
          energy <- energy + penalties[1]
        }
        
        # Checking for consecutive shifts
        if (j < length(shiftstaffing)) # Next shift in same day
        {
          for (k in shiftindices[j]:(shiftindices[j+1]-1))
          {
            if (sum(x[shiftindices[j+1]:(shiftindices[j+2]-1),i]==x[k,i]) > 0)
            {
              energy <- energy + penalties[2]
            }
          }
        }
        else if(i < ndays) # If last shift of the day, following morning
        {
          for (k in shiftindices[j]:(shiftindices[j+1]-1))
          {
            if (sum(x[1:(shiftindices[2]-1),i+1]==x[k,i]) > 0)
            {
              energy <- energy + penalties[2]
            }
          }
        }
        else if(periodic) # If last shift of last day AND periodic, first morning
        {
          for (k in shiftindices[j]:(shiftindices[j+1]-1))
          {
            if (sum(x[1:(shiftindices[2]-1),1]==x[k,i]) > 0)
            {
              energy <- energy + penalties[2]
            }
          }
        }
      }
    }
    shifts <- rep(0, ntotal)
    # Count number of shifts for each worker
    # and check for 3 days in a row
    for (i in 1:ntotal)
    {
      inarow <- 0
      days <- colSums(x==i)
      beginstreakbool <- periodic
      beginstreak <- 0
      for (j in 1:ndays)
      {
        # The streak continues...
        if(days[j] > 0)
        {
          inarow <- inarow + 1
          if (beginstreakbool)
          {
            beginstreak <- beginstreak + 1
          }
          if (inarow == 3)
          {
            energy <- energy + penalties[3]
          }
        }
        else # The streak is over.
        {
          inarow <- 0
          beginstreakbool <- FALSE
        }
      }
      # If periodic and you begin and end the week with partial streaks,
      # check if they connect to give a full streak
      if (periodic && beginstreak < 3 && inarow < 3 && beginstreak + inarow >= 3)
      {
        energy <- energy + penalties[3]
      }
      shifts[i] <- sum(x==i)
    }
    # Penalize energy according to standard deviation among shifts of same worker type
    if (nsupervisors > 1) {
      energy <- energy + penalties[4] * sd(shifts[1:nsupervisors])
    }
    if (nemployees > 1) {
      energy <- energy + penalties[4] * sd(shifts[(nsupervisors+1):ntotal])
    }
    return(energy)
  }
  
  # Retrieve parameters
  nsupervisors <- parameter_df$nsupervisors[index]
  nemployees <- parameter_df$nemployees[index]
  ndays <- parameter_df$ndays[index]
  shiftstaffing <- unlist(parameter_df$shiftstaffing[index])
  penalties <- unlist(parameter_df$penalties[index])
  niter <- parameter_df$niter[index]
  method <- parameter_df$cooling[index]
  step <- parameter_df$step[index]

  # (Some things which are calculated automatically)
  ntotal <- nsupervisors + nemployees
  shiftindices <- c(1,cumsum(shiftstaffing)+1)
  
  # Initial point: put people into the schedule in numerical order each day
  x <- matrix(rep(pmin(1:sum(shiftstaffing),ntotal),ndays), nrow = sum(shiftstaffing), ncol = ndays)
  # Starting energy...
  energy <- rep(0,niter)
  energy[1] <- get_energy(x, penalties)
  
  # Now do the Metropolis loop.
  accepts <- 0
  for (i in 2:niter)
  {
    # Current temperature (from our cooling schedule)
    temperature <- coolingschedule(i, niter, method)
    
    # Choose a random shift and switch a worker
    prop <- x
    
    for (j in 1:step)
    {
      posx <- sample(1:sum(shiftstaffing),1)
      posy <- sample(1:ndays,1)
      prop[posx,posy] <- sample(1:ntotal,1)
      myblock <- min(which(shiftindices>posx))-1
      
      # If the new worker is already on that shift, try again
      while(sum(prop[shiftindices[myblock]:(shiftindices[myblock+1]-1),posy]==prop[posx,posy])>1)
      {
        prop[posx,posy] <- sample(1:ntotal,1)
      }
    }
    
    # Calculate energy of new configuration
    new_energy <- get_energy(prop, penalties)

        # With the appropriate probability, accept and update state/energy
    if(runif(1) < exp((energy[i-1]-new_energy)/temperature))
    {
      energy[i] <- new_energy
      x <- prop
      accepts <- accepts + 1
    }
    else # otherwise, leave things as they stand.
    {
      energy[i] <- energy[i-1]
    }
  }
  # Return the results
  return(c(list(energy), min(energy), (match(min(energy), energy) / niter), tail(energy, n=1), (accepts / (niter - 1))))
}

## END OF FUNCTION SECTION


##-------Experiment and Visualization Section----------------------------------

#Check if greedy approach (constant temperature) works
subset.greedy <- parameter_df %>% filter(cooling == "greedy")
set.seed(553)
subset.greedy <- subset.greedy[sample(nrow(subset.greedy), (dim(subset.greedy)[1] %/% 200)), ]
row.names(subset.greedy) <- NULL

# optional code to time the process
ptm <- proc.time()
for (i in 1:dim(subset.greedy)[1]){
  print(i)
  result <- simuated_annealing(i)
  subset.greedy$energy[i] <- list(result[[1]])
  subset.greedy$min_energy[i] <- result[[2]]
  subset.greedy$first_min_R[i] <- result[[3]]
  subset.greedy$final_energy[i] <- result[[4]]
  subset.greedy$acceptsR[i] <- toString(result[[5]])
}
#show the time
print(proc.time() - ptm)
subset.greedy$sample <- 1:10

# Visualize Energy Change
plot(unlist(subset.greedy$energy[1]), ylim=c(-1,10), 
     main="Change of energy with constant scalar in MCMC",
     xlab="Iterations", ylab="Energy")
points(unlist(subset.greedy$energy[2]), col = "red")
points(unlist(subset.greedy$energy[3]), col = "blue")
legend("topright", title="Parameters",
       c("A","B","C"), fill=c("Black", "Red", "Blue"), horiz=TRUE)

subset.greedy$shiftstaffing <- toString(subset.greedy$shiftstaffing)


# Subsetting parameter_df based on desired conditions
subset.cooling <- parameter_df %>% filter(nsupervisors == 9,
                                         nemployees == 17,
                                         ndays == 14,
                                         niter == 1e+5,
                                         step == 4)

# Get individual results
result.greedy <- simuated_annealing(2700)
result.lin <- simuated_annealing(3240)
result.log <- simuated_annealing(3780)
result.expo <- simuated_annealing(4320)

print(unlist(result.greedy[2:5]))
print(unlist(result.lin[2:5]))
print(unlist(result.log[2:5]))
print(unlist(result.expo[2:5]))

# Plot the energy change overtime for different methods
methods <- ggplot() + 
  geom_line(aes(y = result.greedy[[1]][1:500000], x = 1:500000, colour = "Greedy")) + 
  geom_line(aes(y = result.lin[[1]][1:500000], x = 1:500000, colour = "Linear")) + 
  # geom_line(aes(y = result.log[[1]][1:500000], x = 1:500000, colour = "Logarithmic")) +
  geom_line(aes(y = result.expo[[1]][1:500000], x = 1:500000, colour = "Exponential")) +
  labs(title="Comparison of different methods", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="Methods",
                      values=c(Greedy = "black", Linear="red", Logarithmic="blue", Exponential="green")) +
  theme_dw()
methods  


# Step Size
subset.step <- parameter_df %>% filter(nsupervisors == 9,
                                          nemployees == 17,
                                          ndays == 14,
                                          niter == 1e+5,
                                          cooling == "log")

result.step1 <- simuated_annealing(1339)
result.step4 <- simuated_annealing(3499)
result.step7 <- simuated_annealing(5659)
result.step10 <- simuated_annealing(7819)

# Plot the energy change overtime for different step sizes
stepsizes <- ggplot() + 
  geom_line(aes(y = result.step1[[1]], x = 1:100000, colour = "step_1")) + 
  geom_line(aes(y = result.step4[[1]], x = 1:100000, colour = "step_4")) + 
  geom_line(aes(y = result.step7[[1]], x = 1:100000, colour = "step_7")) +
  geom_line(aes(y = result.step10[[1]], x = 1:100000, colour = "step_10")) +
  labs(title="Comparison of different stepsizes (Simulated Annealing, Logarithmic)", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="stepsizes",
                      values=c(step_1="black", step_4="red", step_7="blue", step_10="green")) +
  theme_classic()
stepsizes  


# Staff vs Schedule
subset.shifts<- parameter_df %>% filter(nsupervisors == 9,
                                          nemployees == 25,
                                          ndays == 7,
                                          niter == 1e+5,
                                          step == 1)

result.5_9log <- simuated_annealing(754)
result.5_9lin <- simuated_annealing(1294)
result.5_9expo <- simuated_annealing(1834)


result.9_25log <- simuated_annealing(761)
result.9_25lin <- simuated_annealing(1301)
result.9_25expo <- simuated_annealing(1841)


small_staff <- ggplot() + 
  geom_line(aes(y = result.5_9log[[1]], x = 1:100000, colour = "Linear")) + 
  geom_line(aes(y = result.5_9lin[[1]], x = 1:100000, colour = "Logarithmic")) + 
  geom_line(aes(y = result.5_9expo[[1]], x = 1:100000, colour = "Exponential")) +
  labs(title="5 supervisors and 9 employees, with 2 people per 5 shifts per day", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="small_staff",
                      values=c(Linear="red", Logarithmic="blue", Exponential="green")) +
  theme_classic()


big_staff <- ggplot() + 
  geom_line(aes(y = result.9_25log[[1]], x = 1:100000, colour = "Linear")) + 
  geom_line(aes(y = result.9_25lin[[1]], x = 1:100000, colour = "Logarithmic")) + 
  geom_line(aes(y = result.9_25expo[[1]], x = 1:100000, colour = "Exponential")) +
  labs(title="9 supervisors and 25 employees, with 2 people per 5 shifts per day", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="big_staff",
                      values=c(Linear="red", Logarithmic="blue", Exponential="green")) +
  theme_classic()

grid.newpage()
grid.draw(rbind(ggplotGrob(small_staff), ggplotGrob(big_staff), size = "last"))


#Constraint and Penalties
subset.constraint<- parameter_df %>% filter(nsupervisors == 5,
                                        nemployees == 17,
                                        ndays == 7,
                                        niter == 1e+5,
                                        step == 1)

result.lin_1 <- simuated_annealing(727)
result.lin_2 <- simuated_annealing(787)
result.lin_3 <- simuated_annealing(847)

result.log_1 <- simuated_annealing(1267)
result.log_2 <- simuated_annealing(1327)
result.log_3 <- simuated_annealing(1387)

result.expo_1 <- simuated_annealing(1807)
result.expo_2 <- simuated_annealing(1867)
result.expo_3 <- simuated_annealing(1927)


lin_pen <- ggplot() + 
  geom_line(aes(y = result.lin_1[[1]], x = 1:100000, colour = "Unbiased"), alpha = 0.7) + 
  geom_line(aes(y = result.lin_2[[1]], x = 1:100000, colour = "Mid_Biased"), alpha = 0.7) + 
  geom_line(aes(y = result.lin_3[[1]], x = 1:100000, colour = "Very_biased"), alpha = 0.7) +
  labs(title="Linear Simulated Annealing on Cost Functions", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="lin_pen",
                      values=c(Unbiased="gold3", Mid_Biased="skyblue3", Very_biased="magenta3")) +
  theme_classic()


log_pen <- ggplot() + 
  geom_line(aes(y = result.log_1[[1]], x = 1:100000, colour = "Unbiased"), alpha = 0.7) + 
  geom_line(aes(y = result.log_2[[1]], x = 1:100000, colour = "Mid_Biased"), alpha = 0.7) + 
  geom_line(aes(y = result.log_3[[1]], x = 1:100000, colour = "Very_biased"), alpha = 0.7) +
  labs(title="Logarithmic Simulated Annealing on Cost Functions", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="log_pen",
                      values=c(Unbiased="gold3", Mid_Biased="skyblue3", Very_biased="magenta3")) +
  theme_classic()

expo_pen <- ggplot() + 
  geom_line(aes(y = result.expo_1[[1]], x = 1:100000, colour = "Unbiased"), alpha = 0.7) + 
  geom_line(aes(y = result.expo_2[[1]], x = 1:100000, colour = "Mid_Biased"), alpha = 0.7) + 
  geom_line(aes(y = result.expo_3[[1]], x = 1:100000, colour = "Very_biased"), alpha = 0.7) +
  labs(title="Exponential Simulated Annealing on Cost Functions", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="expo_pen",
                      values=c(Unbiased="gold3", Mid_Biased="skyblue3", Very_biased="magenta3")) +
  theme_classic()


grid.newpage()
grid.draw(rbind(ggplotGrob(lin_pen), ggplotGrob(log_pen), ggplotGrob(expo_pen), size = "last"))


#Shifts
subset.shifts<- parameter_df %>% filter(nsupervisors == 9,
                                            nemployees == 25,
                                            ndays == 14,
                                            niter == 1e+5,
                                            step == 1,
                                            cooling == "log")


result.shift1 <- simuated_annealing(1282)
result.shift2 <- simuated_annealing(1309)
result.shift3 <- simuated_annealing(1319)

shifts <- ggplot() + 
  geom_line(aes(y = result.shift1[[1]], x = 1:100000, colour = "Shifts_A"), alpha = 0.7) + 
  geom_line(aes(y = result.shift2[[1]], x = 1:100000, colour = "Shifts_B"), alpha = 0.7) + 
  geom_line(aes(y = result.shift3[[1]], x = 1:100000, colour = "Shifts_C"), alpha = 0.7) +
  labs(title="Logarithmic Simulated Annealing on Different Shifts", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="shifts",
                      values=c(Shifts_A="brown4", Shifts_B="deepskyblue1", Shifts_C="olivedrab3")) +
  theme_classic()
shifts


# Number of Iterations
subset.niter <- parameter_df %>% filter(nsupervisors == 5,
                                        nemployees == 25,
                                        ndays == 14,
                                        step == 1)

result.lin_50k <- simuated_annealing(708)
result.lin_100k <- simuated_annealing(888)
result.lin_500k <- simuated_annealing(1068)

result.log_50k <- simuated_annealing(1248)
result.log_100k <- simuated_annealing(1428)
result.log_500k <- simuated_annealing(1608)

result.expo_50k <- simuated_annealing(1788)
result.expo_100k <- simuated_annealing(1968)
result.expo_500k <- simuated_annealing(2148)


lin_iter <- ggplot() + 
  geom_line(aes(y = result.lin_50k[[1]], x = 1:50000, colour = "fifty_thousand")) + 
  geom_line(aes(y = result.lin_100k[[1]], x = 1:100000, colour = "hundred_thousand")) + 
  geom_line(aes(y = result.lin_500k[[1]], x = 1:500000, colour = "five_hundred_thousand")) +
  labs(title="Linear Simulated Annealing with differen max iterations", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="lin_iter",
                      values=c(fifty_thousand ="chartreuse", hundred_thousand="khaki", five_hundred_thousand="black")) +
  theme_classic()


log_iter <- ggplot() + 
  geom_line(aes(y = result.log_50k[[1]], x = 1:50000, colour = "fifty_thousand")) + 
  geom_line(aes(y = result.log_100k[[1]], x = 1:100000, colour = "hundred_thousand")) + 
  geom_line(aes(y = result.log_500k[[1]], x = 1:500000, colour = "five_hundred_thousand")) +
  labs(title="Logarithmic Simulated Annealing with differen max iterations", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="log_iter",
                      values=c(fifty_thousand ="chartreuse", hundred_thousand="khaki", five_hundred_thousand="black")) +
  theme_classic()

expo_iter <- ggplot() + 
  geom_line(aes(y = result.expo_50k[[1]], x = 1:50000, colour = "fifty_thousand")) + 
  geom_line(aes(y = result.expo_100k[[1]], x = 1:100000, colour = "hundred_thousand")) + 
  geom_line(aes(y = result.expo_500k[[1]], x = 1:500000, colour = "five_hundred_thousand")) +
  labs(title="Exponential Simulated Annealing with differen max iterations", x ="Iterations", y = "Energy") +
  scale_colour_manual(name="expo_iter",
                      values=c(fifty_thousand ="chartreuse", hundred_thousand="khaki", five_hundred_thousand="black")) +
  theme_classic()

grid.newpage()
grid.draw(rbind(ggplotGrob(lin_iter), ggplotGrob(log_iter), ggplotGrob(expo_iter), size = "last"))

print(result.lin_50k[[4]])
print(result.lin_100k[[4]])
print(result.lin_500k[[4]])

print(result.log_50k[[4]])
print(result.log_100k[[4]])
print(result.log_500k[[4]])

print(result.expo_50k[[4]])
print(result.expo_100k[[4]])
print(result.expo_500k[[4]])

## END OF EXPERIMENT SECTION --------------------
