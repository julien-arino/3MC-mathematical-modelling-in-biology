# Example simulation of a simple SIS model
library(markovchain)

# Source a file with a few helpful functions for plotting (nice axes labels, crop figure)
source("useful_functions.R")

# Total population
Pop = 5
# Initial number of infectious
I_0 = 2
# Parameters
gamma = 1/5
R_0 = 1.5
# R0 would be (beta/gamma)*S0, so beta=R0*gamma/S0
beta = R_0*gamma/(Pop-I_0)
# Time step
Delta_t = 1

# Make the transition matrix
T = mat.or.vec(nr = (Pop+1), nc = (Pop+1))
for (row in 2:Pop) {
  I = row-1
  mv_right = gamma*I*Delta_t # Recoveries
  mv_left = beta*I*(Pop-I)*Delta_t # Infections
  T[row,(row-1)] = mv_right
  T[row,(row+1)] = mv_left
}
# Last row only has move left
T[(Pop+1),Pop] = gamma*(Pop)*Delta_t
# Check that we don't have too large values
if (max(rowSums(T))>1) {
  T = T/max(rowSums(T))
}
diag(T) = 1-rowSums(T)


mcSIS <- new("markovchain", 
             states = sprintf("I_%d", 0:Pop),
             transitionMatrix = T,
             name = "SIS")

# Show some info about the chain
summary(mcSIS)

# The vector of steady states. Here, all mass should be in I_0 (the first state)
steadyStates(mcSIS)

# Here, we only have one absorbing state, so this is trivial
meanRecurrenceTime(mcSIS)

# Hitting probabilities
hittingProbabilities(mcSIS)

# Mean absorption time
meanAbsorptionTime(mcSIS)
