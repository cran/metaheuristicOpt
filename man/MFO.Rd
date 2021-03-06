% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MFO.Algorithm.R
\name{MFO}
\alias{MFO}
\title{Optimization using Moth Flame Optimizer}
\usage{
MFO(FUN, optimType = "MIN", numVar, numPopulation = 40,
  maxIter = 500, rangeVar)
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
}
\value{
\code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
        and \code{vn} is value of \code{n-th} variable.
}
\description{
This is the internal function that implements Moth Flame Optimization
Algorithm. It is used to solve continuous optimization tasks.
Users do not need to call it directly,
but just use \code{\link{metaOpt}}.
}
\details{
This algorithm was proposed (Mirjalili, 2015). The main inspiration of
this optimizer is the navigation method of moths in nature called transverse
orientation. Moths fly in night by maintaining a fixed angle with respect to
the moon, a very effective mechanism for travelling in a straight line
for long distances. However, these fancy insects are trapped in a useless/deadly
spiral path around artificial lights.

In order to find the optimal solution, the algorithm follow the following steps.
\itemize{
\item Initialization: Initialize the first population of moth randomly,
      calculate the fitness of moth and find the best moth as the best flame obtained so far
      The flame indicate the best position obtained by motion of moth. So in this step, position of
      flame will same with the position of moth.
\item Update Moth Position: All moth move around the corresponding flame.
      In every iteration, the number flame is decreasing over the iteration.
      So at the end of iteration all moth will move around the best solution obtained so far.
\item Replace a flame with the position of moth if a moth becomes fitter than flame
\item Check termination criteria, if termination criterion is satisfied, return the
      best flame as the optimal solution for given problem. Otherwise, back to Update Moth Position steps.
}
}
\examples{
##################################
## Optimizing the schewefel's problem 2.22 function

# define schewefel's problem 2.22 function as objective function
schewefels2.22 <- function(x){
   return(sum(abs(x)+prod(abs(x))))
}

## Define parameter
numVar <- 5
rangeVar <- matrix(c(-10,10), nrow=2)

## calculate the optimum solution using Moth Flame Optimizer
resultMFO <- MFO(schewefels2.22, optimType="MIN", numVar, numPopulation=20,
                 maxIter=100, rangeVar)

## calculate the optimum value using schewefel's problem 2.22 function
optimum.value <- schewefels2.22(resultMFO)

}
\references{
Seyedali Mirjalili, Moth-flame optimization algorithm: A novel nature-inspired
heuristic paradigm, Knowledge-Based Systems, Volume 89, 2015, Pages 228-249,
ISSN 0950-7051, https://doi.org/10.1016/j.knosys.2015.07.006
}
\seealso{
\code{\link{metaOpt}}
}
