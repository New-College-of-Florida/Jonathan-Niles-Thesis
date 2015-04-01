# set the resolution
res <- '200kb'
tf <- 'pol2'

table <- paste(tf, res, 'window.counts.txt', sep='.')
cat("Reading file: ", table, "\n")

# read in ctcf binding from computed bed
ctcf <- read.table(table)

fname <- paste('/home/jniles/thesis/sync/plots/epigenetics/', tf, '.',  res, '.ps', sep='')
cat("Saving to file: ", fname, "\n")

postscript(fname)

# plot ctcf binding in red markers
plot(ctcf[,1], ctcf[,2], col='darkred',xaxt="n", xlab="Distance from Boundary (bp)", ylab="Depth", lty=3)

# adjust labels based on distance to domain boundary
# recall that the window size is 100 base pairs
axis(1, at=seq(0,5000,500), labels=seq(-500000,500000,100000))

# line at domain boundary
abline(v = 2500, col = "gray60", lwd=3, lty=3)

plabel <- paste(toupper(tf), 'Binding', sep=' ')

legend('topleft',
        c(plabel),
        lty=1,
        col=c('darkred'),
        bty='n')

dev.off()
