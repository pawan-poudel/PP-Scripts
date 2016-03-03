\name{scale_colour_Publication}
\alias{scale_colour_Publication}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Additional ggplot settings for colors
}
\description{
Additional ggplot settings for colors
}
\usage{
scale_colour_Publication(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Additional ggplot settings for colors
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (...) 
{
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", 
        "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", 
        "#fb9a99", "#984ea3", "#ffff33")), ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
