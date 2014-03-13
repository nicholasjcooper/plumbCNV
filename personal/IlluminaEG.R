
Illumina.eg <- function (x, mode = "scores") 
{
    valid.mode <- c("date", "data.frame", "md5sum", "scores")
    x <- file(x, open = "rb")
    date <- .read.egt.header(x)
    if (mode == "date") {
        close(x)
        return(date)
    }
    else if (mode %in% valid.mode) {
        cat("Illumina cluster file dated:",date)
        result <- NA
        tryCatch(result <- .read.egt.body(x, read.scores = (mode == 
          "scores")), error = function(e) {
            stop("Couldn't read EGT file", call. = FALSE)
          }, finally = close(x))
        if (mode == "data.frame" || mode == "scores" || is.na(result)) 
            return(result)
        else if (mode == "md5sum") {
            tmp.file <- tempfile()
            .egt.write.text(tmp.file, result)
            md5 <- tools::md5sum(tmp.file)
            unlink(tmp.file)
            return(md5[[1]])
        }
    }
}


.read.egt.body <- function (x, read.scores = FALSE) 
{
    version <- .net.readInt32(x)
    if (version > 7) 
        stop("Unsupported new version:", version, call. = FALSE)
    if (version < 5) 
        stop("Unsupported old version:", version, ". Use 'Broccoli' from outmodedbonsai at http://sourceforge.net/", 
            call. = FALSE)
    .net.readString(x)
    snp.count <- .net.readInt32(x)
    cc.AA <- integer(snp.count)
    cc.AB <- integer(snp.count)
    cc.BB <- integer(snp.count)
    for (i in 1:snp.count) {
        temp <- .read.cluster.definition.record(x, version)
        cc.AA[i] <- temp[1]
        cc.AB[i] <- temp[2]
        cc.BB[i] <- temp[3]
    }
    CSep <- numeric(snp.count)
    total.score <- numeric(snp.count)
    origScore <- numeric(snp.count)
    if (read.scores) {
        for (i in 1:snp.count) {
            CSep[i] <- .net.readFloat(x, 1)
            total.score[i] <- .net.readFloat(x, 1)
            origScore[i] <- .net.readFloat(x, 1)
            .net.readBool(x)
        }
    }
    else for (i in 1:snp.count) {
        .net.readFloat(x, 3)
        .net.readBool(x)
    }
    .net.readString(x, snp.count)
    snp.names <- character(snp.count)
    for (i in 1:snp.count) {
        snp.names[i] <- .net.readString(x)
    }
    if (read.scores) {
        data.frame(row.names = snp.names, SNP = snp.names, AA = cc.AA, 
            AB = cc.AB, BB = cc.BB, Sep = CSep, Score = total.score, 
            origScore = origScore)
    }
    else data.frame(row.names = snp.names, SNP = snp.names, AA = cc.AA, 
        AB = cc.AB, BB = cc.BB)
}


.read.egt.header <- function (x) 
{
    version <- .net.readInt32(x)
    if (version > 3) 
        stop("Invalid cluster file", call. = FALSE)
    gc <- .net.readString(x)
    cluster <- .net.readString(x)
    call <- .net.readString(x)
    normalization <- .net.readString(x)
    create <- .net.readString(x)
    is.WGT <- .net.readBool(x)
    if (version == 2) 
        stop("Not Implemented; Use 'Broccoli' from outmodedbonsai at http://sourceforge.net/", 
            call. = FALSE)
    manifest <- .net.readString(x)
    create
}


.net.readString <- function (con, num = 1L) 
{
    if (num > 1) 
        for (i in 2:num) {
            length <- readBin(con, "integer", n = 1L, size = 1, 
                signed = FALSE)
            if (length > 0) 
                readChar(con, length, useBytes = TRUE)
        }
    length <- readBin(con, "integer", n = 1L, size = 1, signed = FALSE)
    if (length > 0) {
        readChar(con, length, useBytes = TRUE)
    }
    else ""
}


.net.readBool <- function (con, n = 1L) 
{
    readBin(con, "logical", n = n, size = 1, signed = TRUE, endian = "little")
}


.net.readInt32 <- function (con, n = 1L) 
{
    readBin(con, "integer", n = n, size = 4, signed = TRUE, endian = "little")
}


.egt.write.text <- function (file, x) 
{
    write.table(x, file, eol = "\n", col.names = c("SNP", "N_AA", 
        "N_AB", "N_BB"), sep = "\t", quote = FALSE, row.names = FALSE)
}


.net.readFloat <- function (con, n = 1L) 
{
  readBin(con, "numeric", n = n, size = 4, signed = TRUE, endian = "little")
}


.read.cluster.definition.record <- function (x, version) 
{
    result <- .net.readInt32(x, 3)
    if (version >= 6) 
        .net.readFloat(x, 27)
    if (version == 5) 
        .net.readFloat(x, 12)
    result
}


