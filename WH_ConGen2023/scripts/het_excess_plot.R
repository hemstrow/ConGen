library(snpR); library(ggplot2)

eff <- seq(10, 100, by = 10)

res <- data.frame(ne = rep(eff, each = 10), he = 0, ho = 0)
fres <- res

prog <- 1
for(ne in eff){
  cat(ne, "\n")
  pf <- readLines("parmfile")
  nels <- grep("^number_.+_each_population", pf)
  pf[nels] <- paste0("number_", c("females", "males"), "_each_population:\t", ne/2)
  writeLines(pf, "parmfile")
  
  system("/usr/bin/easypop.revised.windows.exe read parmfile", ignore.stdout = TRUE, ignore.stderr = TRUE)
  files <- list.files(pattern = "test.+\\.gen")
  nfiles <- gsub("\\.gen", ".genepop", files)
  file.rename(files, nfiles)
  
  for(file in nfiles){
    cat("\t", file, "\n")
    ts <- read_genepop(file)
    fts <- filter_snps(ts, hwe = 0.001, verbose = FALSE)
    
    ts <- calc_he(ts)
    ts <- calc_ho(ts)
    ts <- get.snpR.stats(ts, stats = c("he", "ho"))$weighted.means
    res[prog,c("he", "ho")] <- unlist(ts[,c("weighted_mean_he", "weighted_mean_ho")])
    
    
    fts <- calc_he(fts)
    fts <- calc_ho(fts)
    fts <- get.snpR.stats(fts, stats = c("he", "ho"))$weighted.means
    fres[prog,c("he", "ho")] <- unlist(fts[,c("weighted_mean_he", "weighted_mean_ho")])
    
    prog <- prog + 1
  }
  file.remove(c(list.files(pattern = "test.+\\.equ"),
                list.files(pattern = "test.+\\.dat"),
                list.files(pattern = "test.+\\.genepop")))
}

tres <- rbind(cbind(res, filtered = FALSE),
              cbind(fres, filtered = TRUE))

ggplot(tres, aes(y = ho-he, x = ne, group = interaction(ne, filtered), color = filtered)) + geom_boxplot() +
  theme_bw() +
  xlab(bquote(N[eb])) +
  ylab("Heterozygote Excess") +
  khroma::scale_color_highcontrast()
  
  


