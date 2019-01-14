# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR <- r*R*(1 - R/K)

    return(list(c(dR)))
  })
}

# The initial state of the system has to be specified as a named vector of state values.
state <- c(R=0.01)

# Parameters has to be specified as a named vector of parameters.
parms <- c(r=1, K=1)

phaseplane(model, state, parms)
