library(photobiologyFlecks)

# Taken on 2021.07.17 in the oat fields of Helsinki university
file.path <- system.file("extdata", "example-oat.csv", package = "photobiologyFlecks")

df <- read.csv(file.path, row.names = 1)

df1 <- df[df$height == "low",]
df2 <- df[df$height == "mid",]
df3 <- df[df$height == "top",]

dfX <- df1
Z <- find_zeros(time = dfX$Time,
                var = dfX$PAR_q)
dfS <- find_flecks(time = dfX$Time,
                   var = dfX$PAR_q,
                   zeroes = Z,
                   minTime = 0, minAmp = 5, minPdiff = 0.05,
                   shadeflecks = FALSE,
                   asmMethod = "max",
                   verbose = FALSE)
plot_ts_fleck(time = dfX$Time,
              var = dfX$PAR_q,
              zeroes = Z,
              fleck.data = dfS)
