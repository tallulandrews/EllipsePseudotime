\name{fit_Ring}
\alias{fit_Ring}
\title{Fits a smoothed ring through 2D data.}
\description{
	Fits an ellipse to 2D data then finds a smoothed path through the highest density of points around the ellipse. Provides a pseudotime for cyclical data starting at a user-defined starting point. 
}
\usage{
 fit_Ring(pcs, start_cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, window_size=20, smoothing.factor=0.7) 
}
\arguments{
  \item{pcs}{ 2D projection of your data, first column is X, second column is Y. }
  \item{start_cell}{ index of the start point for the pseudotime.} 
  \item{direction}{ whether cells progress "clockwise" or "counter" clockwise around the cycle. }
  \item{suppress.plot}{ whether to plot the fitted pseudotime. }
  \item{window_size}{ number of cells to consider when smoothing the distance from the centre of the cycle. }
  \item{smoothing.factor}{ degree of smoothing to apply, higher = more.}
}
\details{
	Fits an ellipse to all data to identify the centre of the cycle. Pseudotime is calculated by walking from the starting cell around the ellipse in the specified direction using a smoothed radius from the centre of the ellipse. Accounts for irregularly shaped cycles.

	Note: smoothing parameters: window_size and smoothing.factor will typically need to be tailored to each dataset.
}
\value{A list of output:
	pseudotime = vector of fit pseudotime values
	fit_x = vector of fit (smoothed) x-coordinates
	fit_y = vector of fit (smoothed) y-coordinates
	pseudo_order = vector of the ordering of points through pseudotime
	ellipse_params = matrix of parameters of the fitted ellipse (in "G" format see conicfit package for details).
	orig_angle = vector of raw angles for each cell around the ellipse.
	rot_angle = vector of angles shifted to go from start cell in the specified direction.
}
\examples{
	# Simulated data:
	set.seed(42)
	params <- matrix(c(2.7, -1.2, 7.16, 4.6, 0.06), ncol=1);
	dat <- calculateEllipse(params, 1000)
	dat[,1] <- dat[,1]+rnorm(1000, sd=0.5)
	dat[,2] <- dat[,2]+rnorm(1000, sd=0.5)
	out <- fit_Ring(dat, suppress.plot=TRUE)
}
\keyword{pseudotime}
