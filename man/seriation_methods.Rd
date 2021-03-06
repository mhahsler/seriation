\name{seriation_methods}
\alias{registry_seriate}
\alias{set_seriation_method}
\alias{get_seriation_method}
\alias{list_seriation_methods}
\alias{show_seriation_methods}
\title{Registry for Seriation Methods}
\description{
A registry to manage methods for seriation.
}
\usage{
list_seriation_methods(kind)
get_seriation_method(kind, name)
set_seriation_method(kind, name, definition, description = NULL, control = list(), ...)

## deprecated
show_seriation_methods(kind)
}
\arguments{
  \item{kind}{the data type the method works on. For example, \code{"dist"},
  \code{"matrix"} or \code{"array"}. If missing, then methods for any type are shown. }
  \item{name}{a short name for the method used to refer to the method in
  \code{seriate()}.}
  \item{definition}{a function containing the method's code.}
  \item{description}{a description of the method. For example, a long name.}
  \item{control}{a list with control arguments and default values.}
  \item{...}{further information that is stored for the method in the
  registry.}
}
\details{
The functions below are convenience function for the registry
\code{registry_seriate}.

\code{list_seriation_method()} lists all available methods for a given
data type (\code{kind}). The result is a vector of character strings
with the short names of the methods. If \code{kind} is missing, then a list of
methods is returned.

\code{get_seriation_method()} returns information (including the
implementing function) about a given method in form of an object of
class \code{"seriation_method"}.

With \code{set_seriation_method()} new seriation methods can be added by the
user. The implementing function (\code{definition}) needs to have the formal
arguments \code{x, control}, where \code{x} is the data object and
\code{control} contains a list with additional information for the method
passed on from \code{seriate()}.  The implementation has to return
a list of objects which can be coerced into \code{ser_permutation_vector}
objects (e.g., integer vectors). The elements in the list
have to be in corresponding order to the dimensions of \code{x}.
}
\author{Michael Hahsler}
\examples{
# Registry
registry_seriate

# List all seriation methods by type
list_seriation_methods()

# List methods for matrix seriation
list_seriation_methods("matrix")

get_seriation_method(name = "BEA")

# Example for defining a new seriation method (reverse identity function for matrix)

# 1. Create the seriation method
seriation_method_reverse <- function(x, control) {
   # return a list of order vectors, one for each dimension
   list(seq(nrow(x), 1), seq(ncol(x), 1))
}

# 2. Register new method
set_seriation_method("matrix", "Reverse", seriation_method_reverse,
    description = "Reverse identity order", control = list())

list_seriation_methods("matrix")
get_seriation_method("matrix", "reverse")

# 3. Use the new seriation methods
seriate(matrix(1:12, ncol=3), "reverse")
}
\keyword{misc}
