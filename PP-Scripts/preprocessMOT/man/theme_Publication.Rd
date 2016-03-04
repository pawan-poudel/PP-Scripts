\name{theme_Publication}
\alias{theme_Publication}
\title{ggplot theme for publications}
\description{ggplot theme for publications}
\usage{
theme_Publication(base_size = 14, base_family = "Helvetica")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base_size}{
font size for the plot
}
  \item{base_family}{
font family
}
}
\details{
Additional setting for the ggplot
}
\references{
Anguraj Sadanandam, et all (2016), MOTassigner
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (base_size = 14, base_family = "Helvetica") 
{
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family) + 
        theme(plot.title = element_text(face = "bold", size = rel(1.2), 
            hjust = 0.5), text = element_text(), panel.background = element_rect(colour = NA), 
            plot.background = element_rect(colour = NA), panel.border = element_rect(colour = NA), 
            axis.title = element_text(face = "bold", size = rel(1)), 
            axis.title.y = element_text(angle = 90, vjust = 2), 
            axis.title.x = element_text(vjust = -0.2), axis.text = element_text(), 
            axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.line = element_line(colour = "black"), axis.ticks = element_line(), 
            panel.grid.major = element_line(colour = "#f0f0f0"), 
            panel.grid.minor = element_blank(), legend.key = element_rect(colour = NA), 
            legend.position = "right", legend.direction = "vertical", 
            legend.key.size = unit(0.5, "cm"), legend.margin = unit(0, 
                "cm"), legend.title = element_text(face = "italic"), 
            plot.margin = unit(c(10, 5, 5, 5), "mm"), strip.background = element_rect(colour = "#f0f0f0", 
                fill = "#f0f0f0"), strip.text = element_text(face = "bold")))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
