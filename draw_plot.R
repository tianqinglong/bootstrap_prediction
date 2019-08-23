library(dplyr)
library(reshape2)
library(ggplot2)

get_the_coverage_prob <- function(output)
{
  bootcp <- sapply(output, function(x) {x$BootCP}) %>% rowMeans(na.rm = T)
  calicp <- sapply(output, function(x) {x$CaliCP}) %>% rowMeans(na.rm = T)
  
  pbcp <- sapply(output, function(x) {x$PBCP}) %>% rowMeans(na.rm = T)
  gpqcp <- sapply(output, function(x) {x$GPQCP}) %>% rowMeans(na.rm = T)
  
  output <- rbind(bootcp, calicp, pbcp, gpqcp) %>% as.data.frame
  names(output) <- c("UpperTail95", "UpperTail90", "LowerTail90", "LowerTail95")
  
  return(output)
}

make_ggplot_data_frame <- function(dataframe, expected_number, pf1)
{
  dm <- dataframe %>% as.matrix %>% melt
  names(dm) <- c("Method", "Quantiles", "Coverage")
  dm$ExpNum <- expected_number
  dm$ProFail <- pf1
  
  return(dm)
}

# ggplot_df %>% ggplot(aes(x = ExpNum, y = Coverage)) + geom_point(aes(col = Method)) + facet_wrap(Quantiles~ProFail)

make_whole_data_frame <- function(Pf1, Er, name)
{
  load(file = name)
  df1 <- get_the_coverage_prob(output)
  df2 <- make_ggplot_data_frame(df1, Er, Pf1)
  return(df2)
}

d1 <- make_whole_data_frame(0.1, 5, "output_05_01.RData")
d2 <- make_whole_data_frame(0.2, 5, "output_05_02.RData")
d3 <- make_whole_data_frame(0.1, 15, "output_15_01.RData")
d4 <- make_whole_data_frame(0.2, 15, "output_15_02.RData")
d5 <- make_whole_data_frame(0.1, 25, "output_25_01.RData")
d6 <- make_whole_data_frame(0.2, 25, "output_25_02.RData")
d7 <- make_whole_data_frame(0.1, 35, "output_35_01.RData")
d8 <- make_whole_data_frame(0.2, 35, "output_35_02.RData")
d9 <- make_whole_data_frame(0.1, 45, "output_45_01.RData")
d10 <- make_whole_data_frame(0.2, 45, "output_45_02.RData")
d11 <- make_whole_data_frame(0.05, 5, "output_05_005.RData")
d12 <- make_whole_data_frame(0.05, 15, "output_15_005.RData")
d13 <- make_whole_data_frame(0.05, 25, "output_25_005.RData")
d14 <- make_whole_data_frame(0.05, 35, "output_35_005.RData")
d15 <- make_whole_data_frame(0.05, 45, "output_45_005.RData")

df <- rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15)
df$hline[df$Quantiles == "UpperTail95"] <- .95
df$hline[df$Quantiles == "LowerTail95"] <- .95
df$hline[df$Quantiles == "UpperTail90"] <- .9
df$hline[df$Quantiles == "LowerTail90"] <- .9

profail_label <- c("Pf1 = 0.05", "Pf1 = 0.1", "Pf1 = 0.2")
names(profail_label) <- c(0.05, 0.1, 0.2)

df %>% ggplot(aes(x = ExpNum, y = Coverage, col = Method)) + geom_point(aes(shape = Method))+geom_line()+facet_grid(ProFail~Quantiles, labeller = labeller(ProFail = profail_label))+
  geom_hline(aes(yintercept = hline), linetype = "dashed") + xlab("Expected Number of Failures")+ylim(c(0.8, 1))+
  ggtitle("Discrete Within Sample Prediction; Pf2 - Pf1 = 0.2; Beta = 2; Random r_star")
