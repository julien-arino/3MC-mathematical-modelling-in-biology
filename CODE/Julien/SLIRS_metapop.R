# SLIRS_metapop
#
# R script that runs integration of the sample p-SLIRS metapopulation
# model in the paper by Arino (2017).
#
# Julien Arino
# February 2017

# Note that without further adaptation, the results provided by this
# script can be expected to be completely unrealistic, the intent is only 
# to provide a workable canvas.


# Load required library. 
# You may need to install this library using install.packages("deSolve").
library(deSolve)


# The vector field
SLIRS_metapop_rhs <- function (t, x, p) {
  # We don't use with(as.list(x)) because we need to name the vectors anyway
  S = x[p$idx_S]
  L = x[p$idx_L]
  I = x[p$idx_I]
  R = x[p$idx_R]
  # We need N as well for the incidence function
  N = S+L+I+R
  # For simplicity, we pre-fill the vector of incidence functions.
  Phi = p$beta*S*I/(N^p$eta)
  # Now we set the values of the derivatives. 
  dS = p$b+p$nu*R-p$d*S-Phi+p$MS%*%S
  dL = Phi-(p$epsilon+p$d)*L+p$ML%*%L
  dI = p$epsilon*L-(p$gamma+p$d)*I+p$MI%*%I
  dR = p$gamma*I-(p$nu+p$d)*R+p$MR%*%R
  dx = list(c(dS, dL, dI, dR))
  return(dx)
}

# Function for computing R0 given we know the 21 entry in V inverse
R0_fct = function(p) {
  Sstar = solve(diag(p$d)-p$MS) %*% p$b
  # Case of mass action incidence
  #F12 = diag(as.vector(p$beta*Sstar))
  # Case of proportional incidence
  F12 = diag(as.vector(p$beta*(Sstar^(1-p$eta))))
  Vinv21 = solve(diag(p$epsilon+p$d)-p$ML) %*% 
    diag(p$epsilon) %*% solve(diag(p$gamma+p$d)-p$MI)
  FVinv = F12 %*% Vinv21
  ev = eigen(FVinv)
  Rzero = max(abs(ev$values))
  return(Rzero)
}

# The following data originates from Bluedot.global and shows the average
# daily number of air passengers flying between each of the countries
# listed below. Average computed over one (unspecified) month.
# pop contains the population of the countries.
pop = c(34.017,1348.932,1224.614,173.593,93.261)*1e6
countries = c("Canada","China","India","Pakistan","Philippines")
T = matrix(data = c(0,1268,900,489,200,
                    1274,0,678,859,150,
                    985,703,0,148,58,
                    515,893,144,0,9,
                    209,174,90,2,0),
           nrow =5, ncol = 5, byrow = TRUE)

# We use a list to store parameters.
p = list()

# In order to compute the movement rates, we use the approximation
# explained in Arino & Portet, JMB 2015.
p$M = mat.or.vec(nr = dim(T)[1], nc = dim(T)[2])
for (from in 1:5) {
  for (to in 1:5) {
    p$M[to,from] = -log(1-T[from,to]/pop[from])
  }
  p$M[from,from] = 0
}
p$MS = p$M-diag(colSums(p$M))
# Assume movement is equal for all compartments.
p$ML = p$MS
p$MI = p$MS
p$MR = p$MS
# Number of patches. (MS should be square.)
p$P = dim(p$MS)[1]

# Set other parameters. Values are chosen mostly arbitrarily.
p$beta = rep(0.1,p$P)
p$nu = rep(0.1,p$P)
p$epsilon = rep(0.01,p$P)
p$d = rep(1/(70*365.25),p$P) # Could be taken based on World Bank data to be more accurate
p$gamma = rep(0.05,p$P)
p$eta = rep(1,p$P)
# We use the numerical trick explained in Arino & Portet, JMB 2015, in
# order to work with a constant total population. Note that this can lead
# to negative birth rates and should only be used if one intends to start
# with an initial condition with the total population equal to its
# equilibrium value.
p$b = (diag(p$d)-p$MS) %*% pop

# Save index of state variable types in state variables vector (we have to use 
# a vector and thus, for instance, the name "S" needs to be defined)
p$idx_S = 1:p$P
p$idx_L = (p$P+1):(2*p$P)
p$idx_I = (2*p$P+1):(3*p$P)
p$idx_R = (3*p$P+1):(4*p$P)
# Set initial conditions. For example, we start with 2 infectious
# individuals in Canada.
L0 = mat.or.vec(p$P,1)
I0 = mat.or.vec(p$P,1)
I0[1] = 2
R0 = mat.or.vec(p$P,1)
S0 = pop-(L0+I0+R0)
# Vector of initial conditions to be passed to ODE solver.
IC = c(S = S0, L = L0, I = I0, R = R0)
# Time span of the simulation (5 years here).
tspan = seq(from = 0, to = 5*365.25, by = 0.1)

# Call of the ODE solver.
sol <- ode(y = IC, times = tspan, func = SLIRS_metapop_rhs, parms = p)

# Put solution in an easier to use form.
times = sol[,"time"]
S = sol[,p$idx_S]
L = sol[,p$idx_L]
I = sol[,p$idx_I]
R = sol[,p$idx_R]
N = S+L+I+R

# A sample (and simple) plot: number infected per 100,000 inhabitants.
xlim = range(times)
ylim = c(0,max(I/N*1e5))
plot(0,
     xlim = xlim, ylim = ylim,
     xlab = "Time (days)", ylab = "Number infectious per 100,000")
for (i in 1:p$P) {
  lines(times,I[,i]/N[,i]*1e5, col = i)
}
legend("topleft",legend = countries,col = 1:p$P, lty = 1)

# Finally, compute the value of R0
R0_fct(p)
