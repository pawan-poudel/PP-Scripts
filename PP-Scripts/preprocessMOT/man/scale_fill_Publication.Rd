\name{scale_fill_Publication}
\alias{scale_fill_Publication}
\title{
Additional ggplot settings for fill
}
\description{
Additional ggplot settings for fill
}
\usage{
scale_fill_Publication(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Setting for the ggolot2
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\references{
Anguraj Sadanandam, et all (2016), MOTassigner
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (...) 
{
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", 
        "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", 
        "#fb9a99", "#984ea3", "#ffff33")), ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
