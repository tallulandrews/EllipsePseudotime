fit_Ring <- function(pcs, start_cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, window_size=20, smoothing.factor=0.7) {

	if (length(start_cell) != 1) {stop("Error: Please supply exactly one starting cell.")}
	if (ncol(pcs) != 2) {warning("Warning: currently only 2D ring fitting is implemented. Using first two columns of pcs.")}

	# Fit Centre
	fit <- conicfit::EllipseDirectFit(pcs)
	fit <- conicfit::AtoG(fit)[[1]]
	fit_centre <- c(fit[1,1], fit[2,1])

	# Calculate Angles
	angles <- atan2(pcs[,2]-fit_centre[2], pcs[,1]-fit_centre[1]);
	if (direction[1] == "clockwise") {
		angles <- pi-angles;
	} else {
		angles <- angles+pi
	}
	converted <- angles-angles[start_cell]
	converted[converted < 0] <- 2*pi+converted[converted <0]

	# Calculate smoothed distance to centre
	ordered_pts <- c(1:length(converted))[order(converted)]
	fit_angles <- converted
	fit_radius <- rep(0, times=length(fit_angles))
	fit_pseudotime <- rep(0, times=length(fit_angles))
	fit_x <- rep(0, times=length(fit_angles))
	fit_y <- rep(0, times=length(fit_angles))
	flank <- ceiling(window_size/2)
	for(i in 1:length(converted)) {
		pos <- ordered_pts[i]
		if (pos-flank < 1) {
			window <- c( (length(ordered_pts)-(flank-pos)):length(ordered_pts),
					 1:pos, (pos+1):(pos+flank))
		} else if (pos+flank > length(ordered_pts)) {
			window <- c( (pos-flank):(pos-1), pos:length(ordered_pts),
					1:( flank-(length(ordered_pts) - pos) ) )
		} else {
			window<- c( (pos-flank):(pos+flank) )
		}

		window_pts <- ordered_pts[window]
		fit_radius[pos] <- median( sqrt((pcs[window_pts,1]-fit_centre[1])^2 + (pcs[window_pts,2]-fit_centre[2])^2  ))
	}

	fit_radius[ordered_pts] <- smooth.spline(as.vector(fit_radius[ordered_pts]), spar=smoothing.factor)$y
	theta <- atan2(pcs[,2]-fit_centre[2], pcs[,1]-fit_centre[1])
	fit_x <- fit_radius*cos(theta) + fit_centre[1]
	fit_y <- fit_radius*sin(theta) + fit_centre[2]

	for(i in 1:length(converted)) {
		pos <- ordered_pts[i]
		if (i == 1) {
			fit_pseudotime[pos] <- 0
		} else {
			pos_prev <- ordered_pts[i-1];
			fit_pseudotime[pos] <- fit_pseudotime[pos_prev]+
				dist(rbind(c(fit_x[pos],fit_y[pos]), 
					c(fit_x[pos_prev],fit_y[pos_prev])));
		}
	}
	
	if (!suppress.plot) {
		vals <- cut(fit_pseudotime, breaks=7)
		cols <- RColorBrewer::brewer.pal(8, "Greys")[2:8]
		plot(pcs[,1], pcs[,2], col=cols[vals], pch=16)
		plot_fit <- calculateEllipse(fit[1,1], fit[2,1], fit[3,1], fit[4,1], fit[5,1], 300)
		lines(plot_fit, lty=3, lwd=1.5, col="grey25")
		points(pcs[start_cell,1], pcs[start_cell,2], pch=1, col="red")
		lines(fit_x[ordered_pts], fit_y[ordered_pts], col="blue", lty=1, lwd=2.5)
	}
	return(list(pseudotime=fit_pseudotime, fit_x=fit_x, fit_y=fit_y, pseudo_order=ordered_pts, ellipse_params=fit, orig_angle=angles, rot_angle=converted))


}

fit_Ring_simple <- function(pcs, start_cell=1, direction=c("clockwise", "counter"), suppress.plot=FALSE, curve_density=360, n_ticks=10, tick_len=0.05) {

	if (length(start_cell) != 1) {stop("Error: Please supply exactly one starting cell.")}
	if (ncol(pcs) != 2) {warning("Warning: currently only 2D ring fitting is implemented. Using first two columns of pcs.")}

	# Fit Centre
	fit <- conicfit::EllipseDirectFit(pcs)
	fit <- conicfit::AtoG(fit)[[1]]
	fit_centre <- c(fit[1,1], fit[2,1])

	# Calculate Angles
	angles <- atan2(pcs[,2]-fit_centre[2], pcs[,1]-fit_centre[1]);
	if (direction[1] == "clockwise") {
		angles <- pi-angles;
	} else {
		angles <- angles+pi
	}
	converted <- angles-angles[start_cell]
	converted[converted < 0] <- 2*pi+converted[converted <0]

	# Calculate smoothed distance to centre
	ordered_pts <- c(1:length(converted))[order(converted)]
	fit_pseudotime <- converted
	
	fit_circle <- conicfit::calculateEllipse(fit[1,1], fit[2,1], fit[3,1], fit[4,1], fit[5,1], curve_density)
	tick_pts <- seq(from=1, to=360, length=n_ticks+1)
	tick_inner <-  conicfit::calculateEllipse(fit[1,1], fit[2,1], fit[3,1]*(1-tick_len), fit[4,1]*(1-tick_len), fit[5,1], curve_density)
	tick_outer <-  conicfit::calculateEllipse(fit[1,1], fit[2,1], fit[3,1]*(1+tick_len), fit[4,1]*(1+tick_len), fit[5,1], curve_density)

	tick_inner <- tick_inner[tick_pts,]
	tick_outer <- tick_outer[tick_pts,]
	colnames(tick_inner) <- c("x", "y")
	colnames(tick_outer) <- c("x", "y")

	if (!suppress.plot) {
		vals <- cut(fit_pseudotime, breaks=7)
		cols <- RColorBrewer::brewer.pal(8, "Greys")[2:8]
		plot(pcs[,1], pcs[,2], col=cols[vals], pch=16)
		lines(fit_circle, lty=3, lwd=1.5, col="grey25")
		points(pcs[start_cell,1], pcs[start_cell,2], pch=1, col="red")
		arrows(tick_inner[,1], tick_inner[,2], tick_outer[,1], tick_outer[,2], length=0, lwd=2)
	}
	return(list(pseudotime=fit_pseudotime, ellipse_pts=fit_circle, pseudo_order=ordered_pts, ellipse_params=fit, ticks_inside=tick_inner, ticks_outside=tick_outer))
}
