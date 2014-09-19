# Example R script to plot scale factor and CC
#
#      by Takanori Nakane (nakane.t@gmail.com)

scale_cc <- read.table("stats.dat", col.names=c("ImageName", "ScaleFactor", "CC"))

# scatter plot
# FIXME: sometimes white lines appear on the plot.
#        seems to be a bug in smoothScatter

smoothScatter(scale_cc$ScaleFactor, scale_cc$CC,
              nbin=200, xlim=c(0, 2), ylim=c(0, 1),
              xlab="Scale Factor", ylab="CC", main="Scale factor & CC")

# histogram

par(mfrow=c(2,2))
plot(scale_cc$ScaleFactor, type='h', ylab="Scale Factor", xlab="#Crystal")
abline(h=mean(scale_cc$ScaleFactor) + (-5:5) * sd(scale_cc$ScaleFactor), col=2)

hist(scale_cc$ScaleFactor, breaks=100, main="Distribution of Scale Factor",
     xlab="Scale Factor")
abline(v=mean(scale_cc$ScaleFactor) + (-5:5) * sd(scale_cc$ScaleFactor), col=2)

plot(scale_cc$CC, type='h', ylab="CC", xlab="#Crystal")
abline(h=mean(scale_cc$CC) + (-5:5) * sd(scale_cc$CC), col=2)

hist(scale_cc$CC, breaks=100,
     main="Distribution of Correlation Coefficient", xlab="CC")
abline(v=mean(scale_cc$CC) + (-5:5) * sd(scale_cc$CC), col=2)
par(mfrow=c(1,1))
