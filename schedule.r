##----------------------------------------------------------------------
## EVERYTHING YOU NEED TO CHANGE IS IN THIS SECTION OR THE FINAL SECTION

# Number of supervisors available
nsupervisors <- 3

# Number of non-supervisor employees
nemployees <- 8

# Number of days in the schedule
ndays <- 7

# TRUE/FALSE: does the schedule repeat? (so for example working
# last night into first morning is penalized)
periodic <- TRUE

# Number of staff present for each shift
# e.g., c(2, 3, 3) has three shifts with 2, 3, 3 people respectively
shiftstaffing <- c(2, 3)

# Penalties for violating each of our four constraints
# Penalty 1: Each shift without a supervisor present
# Penalty 2: Each instance of consecutive shifts for same person (including overnight)
# Penalty 3: Each instance of three days worked in a row
# Penalty 4: Weight for the stdev of each employee's shifts in the schedule
#            (trying to make everyone work similar amounts)
penalties <- c(1, 1, 1, 1)

# Number of iterations (if too small, might not be able to explore enough)
niter <- 250000

# The cooling schedule: what should the temperature be at step i out of n
coolingschedule <- function(i,n)
{
  # If you want greedy optimization, just fix some low temperature
  #return(1e-3)
  
  # Linear cooling schedule from 1 down to 0
  return(1 - (i/(n+1)))
  
  # Logarithmic cooling schedule (goes down really slowly,
  # choose constant to make sure the temp gets low enough)
  #return(0.02/log(1+(i+1)/n))
}

## END OF PARAMETER SECTION
##----------------------------------------------------------------------

# (Some things which are calculated automatically)
ntotal <- nsupervisors + nemployees
shiftindices <- c(1,cumsum(shiftstaffing)+1)

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
  energy <- energy + penalties[4] * (sd(shifts[1:nsupervisors])+sd(shifts[(nsupervisors+1):ntotal]))
  return(energy)
}

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
  temperature <- coolingschedule(i,niter)
  
  # Choose a random shift and switch a worker
  prop <- x
  posx <- sample(1:sum(shiftstaffing),1)
  posy <- sample(1:ndays,1)
  prop[posx,posy] <- sample(1:ntotal,1)
  myblock <- min(which(shiftindices>posx))-1
  
  # If the new worker is already on that shift, try again
  while(sum(prop[shiftindices[myblock]:(shiftindices[myblock+1]-1),posy]==prop[posx,posy])>1)
  {
    prop[posx,posy] <- sample(1:ntotal,1)
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

##----------------------------------------------------------------------
## EVERYTHING YOU NEED TO CHANGE IS IN THIS SECTION OR THE FIRST SECTION

# Various things you might want to look at:

# Plotting the energy over time
plot(energy)

# What was the lowest energy? (Is this necessarily the energy of the final state?)
min(energy)

# Check out our final schedule
print(x)

# Autocorrelation function (not important here since we don't care about
# getting independent samples, we just care about the "final" sample)
acf(energy)

# Acceptance ratio
print(accepts / (niter - 1))

## END OF FINAL SECTION
##----------------------------------------------------------------------