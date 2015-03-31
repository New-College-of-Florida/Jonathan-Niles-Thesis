
# plot ctcf binding
ctcf <- read.table('ctcf.window.counts.txt')
plot(ctcf[,1], ctcf[,2], col='darkred',xaxt = "n",xlab="Distance from Boundary", ylab="Depth")

