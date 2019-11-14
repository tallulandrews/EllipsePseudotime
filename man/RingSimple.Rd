\name{fit_Ring_simple}
\alias{fit_Ring_simple}
\title{Fits angle around an ellispe to 2D data.}
\description{
	Fits an ellipse to 2D data then calculates angles of points around this centre and adjusts them for a user-specified direction from a user-specified start point.
}
\usage{
  fit_Ring_simple(pcs, start_cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, curve_density=360, n_ticks=10, tick_len=0.05) 
}
\arguments{
  \item{pcs}{ 2D projection of your data, first column is X, second column is Y. }
  \item{start_cell}{ index of the start point for the pseudotime.} 
  \item{direction}{ whether cells progress "clockwise" or "counter" clockwise around the cycle. }
  \item{suppress.plot}{ whether to plot the fitted pseudotime. }
  \item{curve_density}{ number of points of the fitted ellipse to calculate for plotting purposes. }
  \item{n_ticks}{ number of ticks to plot around the ring.}
  \item{tick_len}{ length of ticks to plot around the ring.}
}
\details{
	Fits an ellipse to all data to identify the centre of the cycle. Pseudotime is calculated as the angle around the ellipse from a provided starting location in the specified direction.
}
\value{A list of output:
	pseudotime = vector of fit pseudotime values
	ellipse_pts = matrix of x & y coordinates of the fitted ellipse
	pseudo_order = vector of the ordering of points through pseudotime
	ellipse_params = matrix of parameters of the fitted ellipse (in "G" format see conicfit package for details).
	ticks_inside = matix of x & y coordinates of inside ends of the ticks.
	ticks_outside = matix of x & y coordinates of outside ends of the ticks.
}
\examples{
	# Simulated data:
	set.seed(42)
	params <- matrix(c(2.7, -1.2, 7.16, 4.6, 0.06), ncol=1);
	dat <- calculateEllipse(params, 1000)
	dat[,1] <- dat[,1]+rnorm(1000, sd=0.5)
	dat[,2] <- dat[,2]+rnorm(1000, sd=0.5)
	out <- fit_Ring_simple(dat, suppress.plot=TRUE)
}
\keyword{pseudotime}
