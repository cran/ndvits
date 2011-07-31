"AnPhenoMetrics" <-
function (TS, outfile, outgraph= FALSE, Ystart, period, SOSth = 0.5, 
    EOSth = 0.5) 
{
    if (outgraph != FALSE) {
        pdf(outgraph)
    }
    else {
        par(ask = TRUE)
    }
    write(paste("Phenological metrics : \n threshold SOS : ",SOSth, ", EOS : ", EOSth, ", nobsperyear : ", period, sep = ""), outfile, sep = "\t")
    liM = max(TS)
    lim = min(TS)
    final = c()
    ret=c()
    ret$names=rownames(TS)
    colnames(TS)=rep(1:period,length(TS[1,])/period)
    for (j in 1:length(TS[, 1])) {
        res = c()
        par()
        TS[j, ]=as.numeric(TS[j, ])
        plot(ts(as.numeric(TS[j, ]), start = Ystart, freq = period), ylim = c(lim, 
            liM), type = "l", xlab = "time", ylab = "NDVI", main = paste(rownames(TS)[j], "NDVI Time Series"))
        year = c()
        n = 1
        while ((n + period - 1) < length(TS[1, ])) {
            year = rbind(year, TS[j, n:(n + period - 1)])
            n = n + period
        }
        rownames(year) = Ystart:(Ystart + length(year[, 1]) - 
            1)
        meanC = apply(year, 2, mean)
        plot(meanC, ylim = c(lim-0.1, liM), type = "l", xlab = "time", 
            ylab = "NDVI", xaxt = "n", main = paste(rownames(TS)[j], 
                "NDVI Annual Time Series"), lwd = 3)
        axis(side = 1, at = seq(1, period, period/12), labels = month.abb)
        for (i in 1:length(year[, 1])) {
            lines(as.numeric(year[i, ]), col = rainbow(length(year[, 1]))[i])
        }
        legend("bottomright", legend = c("mean", "year"), lwd = c(3, 
            1), col = c("black", "red"))
        meanM = c()
        meanm = c()
        meanC2 = c(meanC[(length(meanC) - 2):length(meanC)], 
            meanC, meanC[1:3])
        for (k in 1:length(meanC)) {
            if (meanC[k] == max(meanC2[(k):(k + 6)]) & meanC[k] > 
                mean(meanC)) 
                meanM = c(meanM, k)
            if (meanC[k] == min(meanC2[(k):(k + 6)]) & meanC[k] < 
                mean(meanC)) 
                meanm = c(meanm, k)
        }
        meanM = sort(meanM)
        meanm = sort(meanm)
        tmp = sort(c(meanM, meanm))
        for (k in 1:(length(tmp) - 1)) {
            if (tmp[k] %in% meanm & tmp[k + 1] %in% meanm) {
                meanm = meanm[-(k - 1 + which.max(c(meanC[tmp[k]], 
                  meanC[tmp[k + 1]])))]
            }
            if (tmp[k] %in% meanM & tmp[k + 1] %in% meanM) 
                meanM = meanM[-(k - 1 + which.min(c(meanC[tmp[k]], 
                  meanC[tmp[k + 1]])))]
        }
        if (length(meanM) != length(meanm)) {
            #cat("Not same number of minimums and maximums, results may be wrong!\n")
            if (length(meanM) > length(meanm)) {
                meanM = meanM[-which.min(meanC[meanM])]
            } else {
                meanm = meanm[-which.max(meanC[meanm])]
            }
        }
        abline(v = c(meanM, meanm), lty = 2)
        text(meanm, -0.1, "min", xpd = NA, cex = 1.2)
        text(meanM, -0.1, "max", xpd = NA, cex = 1.2)
        if (length(meanM) == 1 & length(meanm) == 1 & mean(meanC) > 0.2 & mean(meanC) < 0.7) {
            par(cex=0.9)
            Mdm = meanM
            mdm = meanm
            Mm = max(meanC)
            if (mdm - 2 * (period/12) < 0) {
                st = mdm - 2 * (period/12) + period
                yst = -1
            } else {
                st = mdm - 2 * (period/12)
                yst = 0
            }
            if (mdm + 14 * (period/12) > (2 * period)) {
                ed = mdm + 14 * (period/12) - 2 * period
                yed = 2
            } else {
                ed = mdm + 14 * (period/12) - period
                yed = 1
            }
            for (i in 1:(length(year[, 1]) - 1)) {
                if (yed - yst == 1) {
                  t = c(year[i + yst, st:period], year[i + yed, 
                    1:ed])
                } else {
                  t = c(year[i + yst, st:period], year[i + ((yed + 
                    yst)/2),], year[i + yed, 1:ed])
                }
                #finding key date in the season
                M = max(t[((Mdm - st)%%36 - 2 * (period/12)):(Mdm + 
                  st + 2 * (period/12))])
                Md = which.max(t[((Mdm - st)%%36 - 2 * (period/12)):(Mdm + 
                  st + 2 * (period/12))]) + ((Mdm - st)%%36 - 2 * (period/12)) - 1
                ml = min(t[1:(4 * period/12)])
                mld = which.min(t[1:(4 * period/12)])
                mr = min(t[(12 * period/12):(16 * period/12)]) 
                mrd = which.min(t[(12 * period/12):(16 * period/12)])+ 12 * period/12 -1
                SOSt = ml + SOSth * (M - ml)
                SOS = mld + which(t[mld:Md] > SOSt)[1] - 1
                EOSt = mr + EOSth * (M - mr)
                EOS = Md + which(t[Md:mrd] < EOSt)[1] - 2
                LOS = EOS - SOS
                cumNDVI = sum(t[EOS:SOS])
                cumNDVI2 = sum(t[EOS:SOS]) - (mean(c(t[EOS], t[SOS])) * LOS)
                plot(1:length(t), t, type = "l", 
                  xaxt = "n", ylim = c(0, liM), xlab = "time", 
                  ylab = "SG filtered NDVI", main = paste(rownames(TS)[j], 
                    "NDVI Series", rownames(year)[i], "-", rownames(year)[i + 
                      1]))
                axis(side = 1, at = seq(((st*12)%%period)/12+1, 16*period/12,period/12), labels = rep(month.abb, 2)[((st*12)%/%period+1):((ed*12)%/%period+12)])
                points(x = Md, y = M, col = "green", pch = 3)
                points(x = mld, y = ml, col = "red", pch = 3)
                points(x = mrd, y = mr, col = "red", pch = 3)
                points(x = SOS, y = t[SOS], col = "blue", pch = 3)
                points(x = EOS, y = t[EOS], col = "blue", pch = 3)
                arrows(SOS, 0.05, EOS, 0.05, code = 3, length = 0.1)
                arrows(SOS, -0.5, SOS, t[SOS], length = 0, lty = 2)
                arrows(EOS, -0.5, EOS, t[EOS], length = 0, lty = 2)
                text(SOS, -0.14, "SOS", xpd = NA)
                text(EOS, -0.14, "EOS", xpd = NA)
                text(mean(c(EOS, SOS)), 0.03, paste("cumNDVI =", 
                  as.character(format(cumNDVI, digits = 5)), 
                  sep = ""))
                mld2=(st+mld)%%period
                Md2=(st+Md)%%period
                SOS2=(st+SOS)%%period
                EOS2=(st+EOS)%%period
                res = rbind(res, c(round(ml,3), mld2, round(M,3), Md2, round(SOSt,3), SOS2, 
                  round(EOSt,3), EOS2, LOS, round(cumNDVI, 3), round(cumNDVI2, 3), round(M - Mm, 3)))
            }
            res = as.data.frame(res)
            names(res) = c("\tml", "mld", "M", "Md", "SOSt", "SOS", 
                "EOSt", "EOS", "LOS", "cumNDVI", "cumNDVI-mean", "diffMax")
            rownames(res) = Ystart:(Ystart + length(year[, 1]) - 
                2)
            if (is.null(ret$year)) {
                ret$year=Ystart:(Ystart + length(year[, 1]) - 2)
            }
            ret$mld=cbind(ret$mld,res$mld)
            ret$M=cbind(ret$M,res$M)
            ret$SOS=cbind(ret$SOS,res$SOS)
            ret$EOS=cbind(ret$EOS,res$EOS)
            ret$LOS=cbind(ret$LOS,res$LOS)
            ret$cumNDVI=cbind(ret$cumNDVI,res$cumNDVI)
            write(rownames(TS)[j], outfile, append = TRUE, sep = "")
            write.table(res, outfile, quote = FALSE, row.names = TRUE, append = TRUE, sep = "\t")
        }
        if (length(meanM) == 2 & length(meanm) == 2 & mean(meanC) > 0.2 & mean(meanC) < 
            0.7) {
            par(cex=0.9)
            mdm = meanm[which.min(meanC[meanm])]
            mm = meanC[mdm]
            if (mdm - 2 * (period/12) < 0) {
                st = mdm - 2 * (period/12) + period
                yst = -1
            } else {
                st = mdm - 2 * (period/12)
                yst = 0
            }
            if (mdm + 14 * (period/12) > (2 * period)) {
                ed = mdm + 14 * (period/12) - 2 * period
                yed = 2
            } else {
                ed = mdm + 14 * (period/12) - period
                yed = 1
            }
            if (meanM[1] > mdm) {
                Mdm = meanM[1]
            } else {
                if (meanM[2] > mdm) {
                  Mdm = meanM[2]
                } else {
                  Mdm = meanM[1] + period
                }
            }
            Mdm2 = meanM[!meanM %in% Mdm]
            mdm2 = meanm[!meanm %in% mdm]
            for (i in (1 - yst):(length(year[, 1]) - 1)) {
                if (yed - yst == 1) {
                  t = c(as.numeric(year[i + yst, st:period]), as.numeric(year[i + yed, 
                    1:ed]))
                } else {
                  t = c(as.numeric(year[i + yst, st:period]), as.numeric(year[i + ((yed + 
                    yst)/2),]), as.numeric(year[i + yed, 1:ed]))
                }
                #finding key date in the season
                M1 = max(t[((Mdm - st)%%36 - 2 * (period/12)):((Mdm - st)%%36 
                  + 2 * (period/12))])
                Md1 = which.max(t[((Mdm - st)%%36 - 2 * (period/12)):((Mdm - st)%%36 
                  + 2 * (period/12))]) + ((Mdm - st)%%36 - 2 * (period/12)) - 1
                ml1 = min(t[1 : (4 * (period/12)+1)])
                mld1 = which.min(t[1 : (4 * (period/12+1))])
                M2 = max(t[((Mdm2 - st)%%36 - 2 * (period/12)):((Mdm2 - st)%%36 
                  + 2 * (period/12))])
                Md2 = which.max(t[((Mdm2 - st)%%36 - 2 * (period/12)):((Mdm2 - st)%%36 
                  + 2 * (period/12))]) + ((Mdm2 - st)%%36 - 2 * (period/12)) - 1
                ml2 = min(t[((mdm2 - st)%%36 - 2 * (period/12)):((mdm2 - st)%%36 
                  + 2 * (period/12))])
                mld2 = which.min(t[((mdm2 - st)%%36 - 2 * (period/12)):((mdm2 - st)%%36 
                  + 2 * (period/12))]) + ((mdm2 - st)%%36 - 2 * (period/12)) - 1
		mr = min(t[(12 * (period/12)):(16* (period/12)+1)])
                mrd = which.min(t[(12 * (period/12)):(16* (period/12)+1)])+12 * (period/12) - 1 
                #calculate metrics season 1
                SOS1t = ml1 + SOSth * (M1 - ml1)
                SOS1 = mld1 + which(t[mld1:Md1] > SOS1t)[1] - 1
                EOS1t = ml2 + EOSth * (M1 - ml2)
                EOS1 = Md1 + which(t[Md1:mld2] < EOS1t)[1] - 2
                LOS1 = EOS1 - SOS1
                cumNDVI1 = sum(t[EOS1:SOS1])
                #calculate metrics season 2
                SOS2t = ml2 + SOSth * (M2 - ml2)
                SOS2 = mld2 + which(t[mld2:Md2] > SOS2t)[1] - 1
                EOS2t = mr + EOSth * (M2 - mr)
                EOS2 = Md2 + which(t[Md2:mrd] < EOS2t)[1] - 2
                LOS2 = EOS2 - SOS2
                cumNDVI2 = sum(t[EOS2:SOS2])
                plot(1:length(t), t, type = "l", 
                  xaxt = "n", ylim = c(0, liM), xlab = "time", 
                  ylab = "SG filtered NDVI", main = paste(rownames(TS)[j], 
                    "NDVI Series", rownames(year)[i+yst], "-", rownames(year)[i + yed]))
                if (yed - yst == 1) {
                    axis(side = 1, at = seq(((st*12)%%period)/12+1, 16*period/12,period/12), labels = rep(month.abb, 2)[((st*12)%/%period+1):((ed*12)%/%period+12)])
                } else {
                    axis(side = 1, at = seq(((st*12)%%period)/12+1, 16*period/12,period/12), labels = rep(month.abb, 3)[((st*12)%/%period+1):((ed*12)%/%period+24)])
                }
                points(x = Md1, y = M1, col = "green", pch = 3)
                points(x = mld1, y = ml1, col = "red", pch = 3)
                points(x = Md2, y = M2, col = "green", pch = 3)
                points(x = mld2, y = ml2, col = "red", pch = 3)
                points(x = mrd, y = mr, col = "red", pch = 3)
                points(x = SOS1, y = t[SOS1], col = "blue", pch = 3)
                points(x = EOS1, y = t[EOS1], col = "blue", pch = 3)
                points(x = SOS2, y = t[SOS2], col = "blue", pch = 3)
                points(x = EOS2, y = t[EOS2], col = "blue", pch = 3)
                arrows(SOS1, 0.05, EOS1, 0.05, code = 3, length = 0.1)
                arrows(SOS1, -0.5, SOS1, t[SOS1], length = 0, lty = 2)
                arrows(EOS1, -0.5, EOS1, t[EOS1], length = 0, lty = 2)
                arrows(SOS2, 0.05, EOS2, 0.05, code = 3, length = 0.1)
                arrows(SOS2, -0.5, SOS2, t[SOS2], length = 0, lty = 2)
                arrows(EOS2, -0.5, EOS2, t[EOS2], length = 0, lty = 2)
                text(c(SOS1,SOS2), -0.14, "SOS", xpd = NA)
                text(c(EOS1,EOS2), -0.14, "EOS", xpd = NA)
                text(mean(c(EOS1, SOS1)), 0.03, paste("cumNDVI =", 
                  as.character(format(cumNDVI1, digits = 5)), 
                  sep = ""))
                text(mean(c(EOS2, SOS2)), 0.03, paste("cumNDVI =", 
                  as.character(format(cumNDVI2, digits = 5)), 
                  sep = ""))
                mld12=(st+mld1)%%period
                Md12=(st+Md1)%%period
                mld22=(st+mld2)%%period
                Md22=(st+Md2)%%period
                SOS12=(st+SOS1)%%period
                EOS12=(st+EOS1)%%period
                SOS22=(st+SOS2)%%period
                EOS22=(st+EOS2)%%period
                res = rbind(res, c(round(ml1,3), mld12, round(M1,3), Md12, round(ml2,3), mld22, round(M2,3), Md22, round(SOS1t,3), SOS12, round(EOS1t,3), EOS12, LOS1, round(cumNDVI1, 5), round(SOS2t,3), SOS22, round(EOS2t, 3), EOS22, LOS2, round(cumNDVI2, 5), round(cumNDVI1+cumNDVI2, 5)))
            }
            res = as.data.frame(res)
            names(res) = c("\tml1", "mld1", "M1", "Md1", "ml2", "mld2", "M2", "Md2", "SOS1t", "SOS1", "EOS1t", "EOS1", "LOS1", "cumNDVI1", "SOS2t", "SOS2", "EOS2t", "EOS2", "LOS2", "cumNDVI2", "cumNDVI-tot")
            rownames(res) = rownames(year)[(1 - yst):(length(year[, 1]) - 1)]
            write(rownames(TS)[j], outfile, append = TRUE, sep = "")
            write.table(res, outfile, quote = FALSE, row.names = TRUE, append = TRUE, sep = "\t")
        }
    }
    if (outgraph != FALSE) {
            dev.off()
    }
    return(ret)
}
