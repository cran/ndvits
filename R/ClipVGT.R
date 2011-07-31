"ClipVGT" <-
function (ndvidirectory, region, xlim, ylim, Ystart, Yend, outdirectory, 
    regionout) 
{
    dicperiod = c("01", "11", "21")
    dicmonths = c("01", "02", "03", "04", "05", "06", "07", "08", 
        "09", "10", "11", "12")
    if (!file.exists(substr(outdirectory, start = 1, stop = nchar(outdirectory) - 
        1))) {
        dir.create(outdirectory)
    }
    for (year in seq(Ystart, Yend, 1)) {
        print(paste("Clipping year", as.character(year)))
        for (month in seq(1, 12, 1)) {
            for (period in seq(1, 3, 1)) {
                filein = paste(ndvidirectory, "NDV_", as.character(year), 
                  dicmonths[month], dicperiod[period], "_", region, 
                  "_Extract.tif", sep = "")
                inGrid = readpartGDAL(filein, xlim, ylim)
                fileout = paste(outdirectory, regionout, as.character(year), 
                  "M", dicmonths[month], "P", as.character(period - 
                    1), ".tif", sep = "")
                writeGDAL(inGrid, fileout, drivername = "GTiff", 
                  type = "Byte", mvFlag = 255)
            }
        }
    }
}
