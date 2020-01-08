# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
library(deBif)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    K=q*A
    N0=f*A
    dN <- r*N*(1 - N/K) - E*P*N^2/(N0^2 + N^2)

    return(list(c(dN)))
  })
}

# The initial state of the system has to be specified as a named vector of state values.
state <- c(N=1.0)

# Parameters has to be specified as a named vector of parameters.
parms <- c(A = 0.5, q = 20, E = 0.314, f = 0.474, P = 0.7, r = 0.1)

# phaseplane(model, state, parms)
bifurcation(model, state, parms)
