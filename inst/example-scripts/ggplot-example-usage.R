library(photobiologyFlecks)
library(ggplot2)
set_theme(theme_classic())

# Taken on 2021.07.17 in the oat fields of Helsinki university
file.path <- system.file("extdata", "example-oat.csv", package = "photobiologyFlecks")

df <- read.csv(file.path, row.names = 1)
df$plot.panel <- factor(df$Time %/% 10)

df1 <- df[df$height == "low",]
# df2 <- df[df$height == "mid",]
# df3 <- df[df$height == "top",]

# plot data
ggplot(data = df1) +
  geom_line(aes(Time, PAR_q), linewidth = 0.2, na.rm = TRUE)

Z <- find_zeros(time = df1$Time,
                var = df1$PAR_q)

dfS <- find_flecks(time = df1$Time,
                   var = df1$PAR_q,
                   minTime = 0, minAmp = 5, minPdiff = 0.05,
                   shadeflecks = FALSE,
                   asmMethod = "max",
                   verbose = FALSE)

full.fig <-
  ggplot(data = df1) +
  geom_line(aes(Time, PAR_q), linewidth = 0.2, na.rm = TRUE) +
  geom_point(aes(peakTime, peak), data = dfS, size = 0.5, colour = "coral2", na.rm = TRUE) +
  geom_point(aes(baselineTime1, baseline1), data = dfS, size = 0.5, colour = "chartreuse3", na.rm = TRUE) +
  geom_point(aes(baselineTime2, baseline2), data = dfS, size = 0.5, colour = "chartreuse3", na.rm = TRUE) +
  geom_text(aes(peakTime, peak, label = no), data = dfS, colour = "coral2",
            position = position_nudge(y = 10), size = 2.7, na.rm = TRUE,
            fontface = "bold")

full.fig
full.fig + xlim(0, 5)
full.fig + xlim(5, 10)

full.fig + facet_wrap(facets = vars(plot.panel))
