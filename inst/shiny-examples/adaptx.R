# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
adaptx <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dx = y
    dy = z
    dz = -alpha*z - beta*y - x + x^2

    return(list(c(dx, dy, dz)))
  })
}

# The initial state of the system has to be specified as a named vector of state values.
state <- c(x = 0.0, y = 0.0, z = 0.0)

# Parameters has to be specified as a named vector of parameters.
parms <- c(alpha = 0.5, beta = 1.0)

bifurcation(adaptx, state, parms)

