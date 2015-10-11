\name{criterion_methods}
\alias{registry_criterion}
\alias{set_criterion_method}
\alias{get_criterion_method}
\alias{list_criterion_methods}
\alias{show_criterion_methods}
\title{Registry for Criterion Methods}
\description{
A registry to manage methods to calculate a criterion value given
data and a permutation.
}
\usage{
list_criterion_methods(kind)
show_criterion_methods(kind)
get_criterion_method(kind, name)
set_criterion_method(kind, name, fun, description = NULL, merit = NA, ...)
}
\arguments{
  \item{kind}{the data type the method works on. For example, \code{"dist"},
  \code{"matrix"} or \code{"array"}.}
  \item{name}{a short name for the method used to refer to the method in
  the function \code{criterion()}.}
  \item{fun}{a function containing the method's code.} 
  \item{description}{a description of the method. For example, a long name.} 
  \item{merit}{a boolean indicating if the criterion measure is a merit
  (\code{TRUE}) or a loss (\code{FALSE}) measure.} 
  \item{...}{further information that is stored for the method in the 
  registry.} 
}
\details{
All methods below are convenience methods for the registry named
\code{registry_criterion}.

\code{list_criterion_method()} lists all available methods for a given
data type (\code{kind}). The result is a vector of character strings
with the short names of the methods. 

\code{show_criterion_method()} shows all available methods for a given
data type (\code{kind}) including a description.

\code{get_criterion_method()} returns information (including the
implementing function) about a given method in form of an object of
class \code{"criterion_method"}.

With \code{set_criterion_method()} new criterion methods can be added by the 
user. The implementing function (\code{fun}) needs to have the
formal arguments \code{x, order, ...}, where \code{x} is the data object,
order is an object of class \code{permutation_vector} and \code{...} can
contain additional information for the method passed on from \code{criterion()}.
The implementation has to return the criterion value as a scalar.
}
\author{Michael Hahsler}
\examples{
## the registry
registry_criterion

## use the convenience functions
list_criterion_methods("dist")

show_criterion_methods("dist")

get_criterion_method("dist", "AR_d")

## define a new method

## a function that return sum of the diagonal elements
criterion_method_matrix_foo <- function(x, order, ...) {
if(!is.null(order)) x <- permute(x,order)
    sum(diag(x))
}

## set new method
set_criterion_method("matrix", "foo", criterion_method_matrix_foo, 
    "foo: a useless demo criterion", FALSE)

list_criterion_methods("matrix")

##use all criterion methods (including the new one)
criterion(matrix(1:9, ncol=3))
}
\keyword{misc}
