library("R.utils")

# Plots the TF binding at different windows of resolution
resolutions <- c("100kb", "200kb", "400kb", "800kb", "1000kb")
TF <- 'mafk'
types <- c("peaks", "boundaries")

for (x in 1:length(types)) {
    type <- types[x]

    for (i in 1:length(resolutions)) {
        res <- resolutions[i]

        # get proper file name
        table <- paste(TF, res, type, 'window.counts.txt', sep='.')
        cat("Reading file: ", table, "\n")

        # read in data binding from computed bed
        data <- read.table(table)

        fname <- paste('/home/jniles/thesis/sync/plots/epigenetics/', TF, '.',  res, '.', type, '.ps', sep='')
        cat("Saving to file: ", fname, "\n")

        postscript(fname)

        # plot data binding in red markers
        z <- lowess(data[,1], data[,2], f=1/10)
        plot(data[,1], data[,2], xaxt="n", xlab="Distance from Boundary", ylab="Depth")
        lines(z, type='l', col='darkred', lwd=5)

        title(main=paste(toupper(TF), capitalize(type), sep=" "))

        # adjust labels based on distance to domain boundary
        # recall that the window size is 100 base pairs
        axis(1, at=seq(0,10000,1000), labels=paste(seq(-500,500,100), "kb", sep=""))

        # line at domain boundary
        abline(v = 5000, col = "gray60", lwd=3, lty=3)

        plabel <- paste(toupper(TF), 'Binding', sep=' ')

        legend('topleft',
                c(plabel),
                lty=1,
                col=c('darkred'),
                bty='n')

        dev.off()
    }
}
