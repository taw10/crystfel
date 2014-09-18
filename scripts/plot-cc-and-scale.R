# Example R script to plot scale factor and CC
#
#      by Takanori Nakane (nakane.t@gmail.com)

scale_cc <- read.table("stats.dat", sep=" ", col.names=c("ImageName", "ScaleFactor", "CC"))

# scatter plot
# FIXME: sometimes white lines appear on the plot.
#        seems to be a bug in smoothScatter

smoothScatter(scale_cc$ScaleFactor, scale_cc$CC,
              nbin=200, xlim=c(0, 2), ylim=c(0, 1),
              xlab="Scale Factor", ylab="CC", main="Scale factor & CC")

# histogram

par(mfrow=c(2,1))
plot(scale_cc$ScaleFactor, type='h', ylab="Scale Factor", xlab="#Crystal")
abline(h=mean(sf[,2])+(-5:5)*sd(sf[,2]), col=2)

hist(scale_cc$ScaleFactor, breaks=100, main="Distribution of Scale Factor", xlab="Scale Factor")
abline(v=mean(sf[,2])+(-5:5)*sd(sf[,2]), col=2)
par(mfrow=c(1,1))
