"OverlayTS" <-
function (TS, Ystart, period = 36, title = "NDVI Time Series", 
    outfile = FALSE) 
{
    if (outfile != FALSE) {
        pdf(outfile)
    }
    liM = max(TS)
    lim = min(TS) - 0.1
    if (dim(as.data.frame(TS))[2] > 1) {
        plot(ts(TS[1, ], start = Ystart, freq = period), main = title, 
            ylab = "NDVI", ylim = c(lim, liM), col = rainbow(1)[1])
        for (i in 2:length(TS[, 1])) lines(ts(TS[i, ], start = Ystart, 
            freq = period), col = rainbow(length(TS[, 1]))[i])
        #legend("bottomright", legend = rownames(TS), lwd = 2, 
        #    col = rainbow(length(TS[, 1])))
    }
    else {
        stop("Error : only one time series, no overlay possible.")
    }
    if (outfile != FALSE) {
        dev.off()
    }
}
"PlotTS" <-
function (TS, outfile = FALSE, Ystart, period, title = "NDVI time series", 
    nb = NULL )
{
    liM = max(TS)+0.05
    lim = min(TS)
    if (outfile != FALSE) {
        pdf(outfile)
    }
    if (is.null(nb)) {
    nb=ifelse(length(TS[, 1])>7, 5, length(TS[, 1]))
    }
    if (length(TS[, 1]) > nb & outfile == FALSE) {
        par(ask = TRUE)
    }
    par(mfrow = c(nb, 1), mar = c(0, 4.1, 0, 
        2.1), oma = c(8, 2, 8, 2))
    for (i in 1:length(TS[, 1])) {
        plot(ts(as.numeric(TS[i, ]), start = Ystart, freq = period), 
            ylim = c(lim, liM), type = "l", xlab = "", xaxt = "n", 
            , yaxt = "n", ylab = "")
        axis(side = ifelse(i%%2 == 1, 2, 4), at = seq(round(lim, 
            1), round(liM, 1), 0.1), labels = seq(round(lim, 
            1), round(liM, 1), 0.1))
        mtext("NDVI", side = ifelse(i%%2 == 1, 2, 4), line = 3, 
            cex = 0.7)
        text((Ystart + length(TS[1, ])/(2 * period)), liM - 
            0.025, rownames(TS)[i], xpd = "NA", cex = 1)
        if (i%%nb == 0 || i == length(TS[, 1])) {
            axis(side = 1, at = seq(Ystart, 2009, 1), labels = seq(Ystart, 
                2009, 1))
            mtext(title, side = 3, line = 3, outer = TRUE)
        }
    }
    if (outfile != FALSE) {
        dev.off()
    }
}
