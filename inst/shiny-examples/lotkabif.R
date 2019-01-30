# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N

    return(list(c(dR, dN)))
  })
}

# The initial state of the system has to be specified as a named vector of state values.
state <- c(R=1, N=0.01)

# Parameters has to be specified as a named vector of parameters.
parms <- c(r=1, K=1, a=1, c=1, delta=0.5)

bifurcation(model, state, parms)

