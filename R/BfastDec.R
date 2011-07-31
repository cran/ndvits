"BfastDec" <-
function (TS, Ystart, period, param = c(), outfile = FALSE, ...) 
{
    if (!require(bfast)) {
        return()
    }
    if (is.null(param$season)) 
        param$season = "harmonic"
    if (is.null(param$h)) 
        param$h = period * 2/(length(as.numeric(TS[1, ])))
    if (is.null(param$maxi)) 
        param$maxi = 10
    if (outfile != FALSE) {
        pdf(outfile)
    }
    else {
        par(ask = TRUE)
    }
    seas = list()
    trend = list()
    for (i in 1:length(TS[, 1])) {
        fit = bfast(ts(as.numeric(TS[i, ]), start = Ystart, freq = period), 
            h = param$h, season = param$season, max.iter = param$maxi, 
            ...)
        plot(fit)
        seas = c(seas, Ystart + as.numeric(fit$output[[1]]$bp.Wt$breakpoints)/period)
        trend = c(trend, Ystart + as.numeric(fit$output[[1]]$bp.Vt$breakpoints)/period)
    }
    if (outfile != FALSE) {
        dev.off()
    }
    res = c()
    res$seas = seas
    res$trend = trend
    return(res)
}
