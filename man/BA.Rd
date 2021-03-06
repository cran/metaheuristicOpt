% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BA.Algorithm.R
\name{BA}
\alias{BA}
\title{Optimization using Bat Algorithm}
\usage{
BA(FUN, optimType = "MIN", numVar, numPopulation = 40, maxIter = 500,
  rangeVar, maxFrequency = 0.1, minFrequency = -0.1, gama = 1,
  alphaBA = 0.1)
}
\arguments{
\item{FUN}{an objective function or cost function,}

\item{optimType}{a string value that represent the type of optimization.
There are two option for this arguments: \code{"MIN"} and \code{"MAX"}.
The default value is \code{"MIN"}, which the function will do minimization.
Otherwise, you can use \code{"MAX"} for maximization problem.
The default value is \code{"MIN"}.}

\item{numVar}{a positive integer to determine the number variables.}

\item{numPopulation}{a positive integer to determine the number populations. The default value is 40.}

\item{maxIter}{a positive integer to determine the maximum number of iterations. The default value is 500.}

\item{rangeVar}{a matrix (\eqn{2 \times n}) containing the range of variables,
where \eqn{n} is the number of variables, and first and second rows
are the lower bound (minimum) and upper bound (maximum) values, respectively.
If all variable have equal upper bound, you can define \code{rangeVar} as
matrix (\eqn{2 \times 1}).}

\item{maxFrequency}{a numeric to determine maximum frequency. The default value is 0.1.}

\item{minFrequency}{a numeric to determine minimum frequency. The default value is -0.1.}

\item{gama}{a numeric greater than equal to 1. It use to increase pulse rate. The default value is 1.}

\item{alphaBA}{a numeric between 0 and 1. It use to decrease loudness. The default value is 0.1.}
}
\value{
\code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
        and \code{vn} is value of \code{n-th} variable.
}
\description{
This is the internal function that implements Bat
Algorithm. It is used to solve continuous optimization tasks.
Users do not need to call it directly,
but just use \code{\link{metaOpt}}.
}
\details{
This algorithm was proposed by (Yang, 2011). It was inspired
by echolocation of bats. Candidate solutions in bat algorithm
are represented by bat. They have flying speed, pulse rate,
loudness and pulse frequency and they move based on them.

In order to find the optimal solution, the algorithm follow the following steps.
\itemize{
\item initialize population randomly.
\item move every candidate solutions based on velocity and pulse frequnecy.
\item move some candidate solutions near globak best randomly.
\item If a candidate solution have better fitness than global best replace it with new
      random candidate solution then increase its pulse rate and decrease its loudness.
\item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
      exit the loop, else back to move every candidate solutions.
}
}
\examples{
##################################
## Optimizing the schewefel's problem 1.2 function

# define schewefel's problem 1.2 function as objective function
schewefels1.2 <- function(x){
  dim <- length(x)
  result <- 0
    for(i in 1:dim){
       result <- result + sum(x[1:i])^2
   }
  return(result)
}

## Define parameter
numVar <- 5
rangeVar <- matrix(c(-10,10), nrow=2)

## calculate the optimum solution using bat algorithm
resultBA <- BA(schewefels1.2, optimType="MIN", numVar, numPopulation=20,
                 maxIter=100, rangeVar)

## calculate the optimum value using schewefel's problem 1.2 function
optimum.value <- schewefels1.2(resultBA)

}
\references{
Yang, X. S., (2011), Bat Algorithm for Multiobjective Optimization,
Int. J. Bio-Inspired Computation, Vol. 3, No. 5, pp.267-274.
}
\seealso{
\code{\link{metaOpt}}
}
