filenames <- list.files(path = ".", pattern = "^out.*\\.csv")

par(ask = TRUE)
for (outfile in filenames) {
    data <- read.csv(outfile)[, -1]
    matplot(data[, 1], log10(data[, -1]), type = 'l', main = outfile,
            xlab = 't\'', ylab = 'B')
}
