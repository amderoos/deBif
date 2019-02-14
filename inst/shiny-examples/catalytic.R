# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model

# Catalytic oscillator
catalytic <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    z <- (1 - x - y - s)
    dx <- 2*q1*(z^2) - 2*q5*(x^2) - q3*x*y
    dy <- q2*z - q6*y - q3*x*y
    ds <- q4*z - k*q4*s

    return(list(c(dx, dy, ds)))
  })
}

# The initial state of the system has to be specified as a named vector of state values.
state <- c(x = 0.0029538, y = 0.76211, s = 0.16781)

# Parameters has to be specified as a named vector of parameters.
parms <- c(q1 = 2.5, q2 = 1.4707, q3 = 10, q4 = 0.0675, q5 = 1, q6 = 0.1, k = 0.4)

bifurcation(catalytic, state, parms)

