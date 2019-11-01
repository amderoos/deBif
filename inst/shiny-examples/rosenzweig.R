# The initial state of the system has to be specified as a named vector of state values.
state <- c(R = 0.05, C = 0.1)

# Parameters has to be specified as a named vector of parameters.
parms <- c(r = 0.5, K = 0.1, a = 5.0, h = 3.0, eps = 0.5, mu = 0.05)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
rosenzweig <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR = r*R*(1 - R/K) - a*R*C/(1 + a*h*R)
    dC = eps*a*R*C/(1 + a*h*R) - mu*C

    return(list(c(dR, dC)))
  })
}

bifurcation(rosenzweig, state, parms)

