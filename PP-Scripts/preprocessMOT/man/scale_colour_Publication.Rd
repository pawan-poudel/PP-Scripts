\name{scale_colour_Publication}
\alias{scale_colour_Publication}
\title{
Additional ggplot settings for colors
}
\description{
Additional ggplot settings for colors
}
\usage{
scale_colour_Publication(...)
}
\arguments{
  \item{\dots}{
   extra parameters that you want to put in the function
}
}
\details{
Additional ggplot settings for colors
}
\references{Anguraj Sadanandam, et all (2016), MOTassigner}
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
