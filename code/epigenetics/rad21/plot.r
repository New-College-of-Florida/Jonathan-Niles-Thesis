

# read in ctcf binding from computed bed
ctcf <- read.table('rad21.window.counts.txt')

postscript('/home/jniles/thesis/sync/plots/epigenetics/rad21.ps')

# plot ctcf binding in red markers
plot(ctcf[,1], ctcf[,2], col='darkred',xaxt = "n",xlab="Distance from Boundary (bp)", ylab="Depth")

# adjust labels based on distance to domain boundary
# recall that the window size is 100 base pairs
axis(1, at=seq(0,5000,500), labels=seq(-500000,500000,100000))

# line at domain boundary
abline(v = 2500, col = "gray60", lwd=3, lty=3)

legend('topleft',
        c("CTCF Binding"),
        lty=1,
        col=c('darkred'),
        bty='n')

dev.off()
