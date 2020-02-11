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
show_seriation_methods(kind)
get_seriation_method(kind, name)
set_seriation_method(kind, name, definition, description = NULL, control = list(), ...)
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
with the short names of the methods.

\code{show_seriation_method()} shows all available methods including
a description.

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
## registry
registry_seriate

## convenience functions
show_seriation_methods()
show_seriation_methods("matrix")

list_seriation_methods("matrix")

get_seriation_method(name = "BEA")

## define a new method

## create a identity function which returns the identity order
seriation_method_identity <- function(x, control) {
   lapply(dim(x), seq)
}

## set new method
set_seriation_method("matrix", "identity", seriation_method_identity,
    "Identity order")

set_seriation_method("array", "identity", seriation_method_identity,
    "Identity order")

show_seriation_methods("matrix")

##use all criterion methods (including the new one)
seriate(matrix(1:12, ncol=3), "identity")
}
\keyword{misc}
