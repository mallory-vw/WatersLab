#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_FragLength_plot.RData")

#load packages
library("ggplot2")
library("reshape2")
library("plyr")
library("scales")

#####load data#####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/Sir2A-KO_Mapping/FragLength/")

ATAC_507_4h <- read.table("ATAC_507_4h_exp1_TLEN.txt")
  colnames(ATAC_507_4h) <- c("TLEN")
  ATAC_507_4h_counts <- count(abs(ATAC_507_4h))
  ATAC_507_4h_counts$logfreq <- log(ATAC_507_4h_counts$freq)
  ATAC_507_4h_counts$proportionfreq <- ATAC_507_4h_counts$freq/(sum(ATAC_507_4h_counts$freq))

ATAC_507_10h <- read.table("ATAC_507_10h_exp1_TLEN.txt")
  colnames(ATAC_507_10h) <- c("TLEN")
  ATAC_507_10h_counts <- count(abs(ATAC_507_10h))
  ATAC_507_10h_counts$logfreq <- log(ATAC_507_10h_counts$freq)
  ATAC_507_10h_counts$proportionfreq <- ATAC_507_10h_counts$freq/(sum(ATAC_507_10h_counts$freq))  
  
ATAC_507_16h <- read.table("ATAC_507_16h_exp1_TLEN.txt")
  colnames(ATAC_507_16h) <- c("TLEN")
  ATAC_507_16h_counts <- count(abs(ATAC_507_16h))
  ATAC_507_16h_counts$logfreq <- log(ATAC_507_16h_counts$freq)
  ATAC_507_16h_counts$proportionfreq <- ATAC_507_16h_counts$freq/(sum(ATAC_507_16h_counts$freq))
  
ATAC_507_22h <- read.table("ATAC_507_22h_exp1_TLEN.txt")
  colnames(ATAC_507_22h) <- c("TLEN")
  ATAC_507_22h_counts <- count(abs(ATAC_507_22h))
  ATAC_507_22h_counts$logfreq <- log(ATAC_507_22h_counts$freq)
  ATAC_507_22h_counts$proportionfreq <- ATAC_507_22h_counts$freq/(sum(ATAC_507_22h_counts$freq))
  
ATAC_507_gDNA <- read.table("ATAC_507_gDNA_exp1_TLEN.txt")
  colnames(ATAC_507_gDNA) <- c("TLEN")
  ATAC_507_gDNA_counts <- count(abs(ATAC_507_gDNA))
  ATAC_507_gDNA_counts$logfreq <- log(ATAC_507_gDNA_counts$freq)
  ATAC_507_gDNA_counts$proportionfreq <- ATAC_507_gDNA_counts$freq/(sum(ATAC_507_gDNA_counts$freq))
  
ATAC_1022_4h <- read.table("ATAC_1022_4h_exp1_TLEN.txt")
  colnames(ATAC_1022_4h) <- c("TLEN")
  ATAC_1022_4h_counts <- count(abs(ATAC_1022_4h))
  ATAC_1022_4h_counts$logfreq <- log(ATAC_1022_4h_counts$freq)
  ATAC_1022_4h_counts$proportionfreq <- ATAC_1022_4h_counts$freq/(sum(ATAC_1022_4h_counts$freq))
  
ATAC_1022_10h <- read.table("ATAC_1022_10h_exp1_TLEN.txt")
  colnames(ATAC_1022_10h) <- c("TLEN")
  ATAC_1022_10h_counts <- count(abs(ATAC_1022_10h))
  ATAC_1022_10h_counts$logfreq <- log(ATAC_1022_10h_counts$freq)
  ATAC_1022_10h_counts$proportionfreq <- ATAC_1022_10h_counts$freq/(sum(ATAC_1022_10h_counts$freq))  
  
ATAC_1022_16h <- read.table("ATAC_1022_16h_exp1_TLEN.txt")
  colnames(ATAC_1022_16h) <- c("TLEN")
  ATAC_1022_16h_counts <- count(abs(ATAC_1022_16h))
  ATAC_1022_16h_counts$logfreq <- log(ATAC_1022_16h_counts$freq)
  ATAC_1022_16h_counts$proportionfreq <- ATAC_1022_16h_counts$freq/(sum(ATAC_1022_16h_counts$freq))
  
ATAC_1022_22h <- read.table("ATAC_1022_22h_exp1_TLEN.txt")
  colnames(ATAC_1022_22h) <- c("TLEN")
  ATAC_1022_22h_counts <- count(abs(ATAC_1022_22h))
  ATAC_1022_22h_counts$logfreq <- log(ATAC_1022_22h_counts$freq)
  ATAC_1022_22h_counts$proportionfreq <- ATAC_1022_22h_counts$freq/(sum(ATAC_1022_22h_counts$freq))
  
ATAC_1022_gDNA <- read.table("ATAC_1022_gDNA_exp1_TLEN.txt")
  colnames(ATAC_1022_gDNA) <- c("TLEN")
  ATAC_1022_gDNA_counts <- count(abs(ATAC_1022_gDNA))
  ATAC_1022_gDNA_counts$logfreq <- log(ATAC_1022_gDNA_counts$freq)
  ATAC_1022_gDNA_counts$proportionfreq <- ATAC_1022_gDNA_counts$freq/(sum(ATAC_1022_gDNA_counts$freq)) 
  
  
#####merge data into 1 data frame#####
mylist <- list(ATAC_507_4h = ATAC_507_4h_counts, 
               ATAC_507_10h = ATAC_507_10h_counts,
               ATAC_507_16h = ATAC_507_16h_counts,
               ATAC_507_22h = ATAC_507_22h_counts,
               ATAC_507_gDNA = ATAC_507_gDNA_counts,
               ATAC_1022_4h = ATAC_1022_4h_counts, 
               ATAC_1022_10h = ATAC_1022_10h_counts,
               ATAC_1022_16h = ATAC_1022_16h_counts,
               ATAC_1022_22h = ATAC_1022_22h_counts,
               ATAC_1022_gDNA = ATAC_1022_gDNA_counts)
  
all_counts <- do.call("rbind", mylist)
all_counts$ID <- rep(names(mylist), sapply(mylist, nrow))
all_counts$ID <- as.factor(all_counts$ID)

mylist_gDNA <- list(ATAC_507_gDNA = ATAC_507_gDNA_counts,
                    ATAC_1022_gDNA = ATAC_1022_gDNA_counts)

gDNA_counts <- do.call("rbind", mylist_gDNA)
gDNA_counts$ID <- rep(names(mylist_gDNA), sapply(mylist_gDNA, nrow))
gDNA_counts$ID <- as.factor(gDNA_counts$ID)

mylist_ATAC <- list(ATAC_507_4h = ATAC_507_4h_counts, 
                    ATAC_507_10h = ATAC_507_10h_counts,
                    ATAC_507_16h = ATAC_507_16h_counts,
                    ATAC_507_22h = ATAC_507_22h_counts,
                    ATAC_1022_4h = ATAC_1022_4h_counts, 
                    ATAC_1022_10h = ATAC_1022_10h_counts,
                    ATAC_1022_16h = ATAC_1022_16h_counts,
                    ATAC_1022_22h = ATAC_1022_22h_counts)


ATAC_counts <- do.call("rbind", mylist_ATAC)
ATAC_counts$ID <- rep(names(mylist_ATAC), sapply(mylist_ATAC, nrow))
ATAC_counts$ID <- as.factor(ATAC_counts$ID)



#####plots#####
all_plot <- ggplot(data = all_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
    geom_line(size = 0.5) +
    scale_colour_brewer(palette="Spectral") +
    ylab("Proportion of Reads \n(frament length < 500bp)") + 
    xlab("Insert Size (bp)") +
    scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
      theme(axis.title.y = element_text(angle = 90, size=20)) +
      theme(axis.title.x = element_text(angle = 0, size=20)) +
      theme(axis.text.y = element_text(size=16)) +
      theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
      theme(legend.title = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))


gDNA_plot <- ggplot(data = gDNA_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
  geom_line(size = 0.5) +
  # scale_colour_brewer(palette="Spectral") +
  ylab("Proportion of Reads \n(frament length < 500bp)") + 
  xlab("Insert Size (bp)") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ATAC_plot <- ggplot(data = ATAC_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
  geom_line(size = 0.5) +
  scale_colour_brewer(palette="Spectral") +
  ylab("Proportion of Reads \n(frament length < 500bp)") + 
  xlab("Insert Size (bp)") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


#####testing sebastians previous plot#####
#####load data####
setwd("J:/III/Waters/Group Members/Mallory/PbergheiUpstreamRegionsSebastian/Testing/FragLengthPlotTest/")

Test_50_30min_S3_R1 <- read.table("50-30min_S3_R1_001_TLEN.txt")
  colnames(Test_50_30min_S3_R1) <- c("TLEN")
  Test_50_30min_S3_R1_counts <- count(abs(Test_50_30min_S3_R1))
  Test_50_30min_S3_R1_counts$logfreq <- log(Test_50_30min_S3_R1_counts$freq)
  Test_50_30min_S3_R1_counts$proportionfreq <- Test_50_30min_S3_R1_counts$freq/(sum(Test_50_30min_S3_R1_counts$freq))

Test_500_30min_S16_R1 <- read.table("500-30min_S16_R1_001_TLEN.txt")
  colnames(Test_500_30min_S16_R1) <- c("TLEN")
  Test_500_30min_S16_R1_counts <- count(abs(Test_500_30min_S16_R1))
  Test_500_30min_S16_R1_counts$logfreq <- log(Test_500_30min_S16_R1_counts$freq)
  Test_500_30min_S16_R1_counts$proportionfreq <- Test_500_30min_S16_R1_counts$freq/(sum(Test_500_30min_S16_R1_counts$freq))
  
Test_5000_30min_S16_R1 <- read.table("5000-30min_S16_R1_001_TLEN.txt")
  colnames(Test_5000_30min_S16_R1) <- c("TLEN")
  Test_5000_30min_S16_R1_counts <- count(abs(Test_5000_30min_S16_R1))
  Test_5000_30min_S16_R1_counts$logfreq <- log(Test_5000_30min_S16_R1_counts$freq)
  Test_5000_30min_S16_R1_counts$proportionfreq <- Test_5000_30min_S16_R1_counts$freq/(sum(Test_5000_30min_S16_R1_counts$freq))
  
Test_G1142_I2hpi_S13_R1 <- read.table("G1142-I2hpi_S13_R1_001_TLEN.txt")
  colnames(Test_G1142_I2hpi_S13_R1) <- c("TLEN")
  Test_G1142_I2hpi_S13_R1_counts <- count(abs(Test_G1142_I2hpi_S13_R1))
  Test_G1142_I2hpi_S13_R1_counts$logfreq <- log(Test_G1142_I2hpi_S13_R1_counts$freq)
  Test_G1142_I2hpi_S13_R1_counts$proportionfreq <- Test_G1142_I2hpi_S13_R1_counts$freq/(sum(Test_G1142_I2hpi_S13_R1_counts$freq))

  
#####merge data into 1 data frame#####
Test_mylist <- list(Test_50_30min_S3_R1 = Test_50_30min_S3_R1_counts, 
                    Test_500_30min_S16_R1 = Test_500_30min_S16_R1_counts,
                    Test_5000_30min_S16_R1 = Test_5000_30min_S16_R1_counts,
                    Test_G1142_I2hpi_S13_R1 = Test_G1142_I2hpi_S13_R1_counts)
  
Test_all_counts <- do.call("rbind", Test_mylist)
  Test_all_counts$ID <- rep(names(Test_mylist), sapply(Test_mylist, nrow))
  Test_all_counts$ID <- as.factor(Test_all_counts$ID)
  
Test_mylist_gDNA <- list(Test_G1142_I2hpi_S13_R1 = Test_G1142_I2hpi_S13_R1_counts)
  
Test_gDNA_counts <- do.call("rbind", Test_mylist_gDNA)
  Test_gDNA_counts$ID <- rep(names(Test_mylist_gDNA), sapply(Test_mylist_gDNA, nrow))
  Test_gDNA_counts$ID <- as.factor(Test_gDNA_counts$ID)
  
Test_mylist_ATAC <- list(Test_50_30min_S3_R1 = Test_50_30min_S3_R1_counts, 
                      Test_500_30min_S16_R1 = Test_500_30min_S16_R1_counts,
                      Test_5000_30min_S16_R1 = Test_5000_30min_S16_R1_counts)
  
  Test_ATAC_counts <- do.call("rbind", Test_mylist_ATAC)
  Test_ATAC_counts$ID <- rep(names(Test_mylist_ATAC), sapply(Test_mylist_ATAC, nrow))
  Test_ATAC_counts$ID <- as.factor(Test_ATAC_counts$ID)  

#####plots#####
Test_all_plot <- ggplot(data = Test_all_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
  geom_line(size = 0.5) +
  scale_colour_brewer(palette="Spectral") +
  ylab("Proportion of Reads \n(frament length < 500bp)") + 
  xlab("Insert Size (bp)") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


Test_gDNA_plot <- ggplot(data = Test_gDNA_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
  geom_line(size = 0.5) +
  # scale_colour_brewer(palette="Spectral") +
  ylab("Proportion of Reads \n(frament length < 500bp)") + 
  xlab("Insert Size (bp)") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
Test_ATAC_plot <- ggplot(data = Test_ATAC_counts, aes(x = TLEN, y = proportionfreq, colour = ID)) +
  geom_line(size = 0.5) +
  # scale_colour_brewer(palette="Spectral") +
  ylab("Proportion of Reads \n(frament length < 500bp)") + 
  xlab("Insert Size (bp)") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(1,500), breaks = c(1,50,100,150,200,250,300,350,400,450,500)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, size=20)) +
  theme(axis.title.x = element_text(angle = 0, size=20)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0)) +
  theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

#####plots to compare#####
all_plot
gDNA_plot
ATAC_plot

Test_all_plot
Test_gDNA_plot
Test_ATAC_plot

#####data metrics#####
#do this on Test_50_30min_S3_R1_counts, compare to Galaxy metrics
#TLEN   freq   logfreq proportionfreq

Test_50_30min_S3_R1_counts_metric <- subset(Test_50_30min_S3_R1_counts, Test_50_30min_S3_R1_counts$TLEN >= 20 & Test_50_30min_S3_R1_counts$TLEN < 1930942)
median(Test_50_30min_S3_R1_counts_metric$TLEN)
sd(Test_50_30min_S3_R1_counts_metric$TLEN)
mean(Test_50_30min_S3_R1_counts_metric$TLEN)


#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Pberghei_upstream_region_FragLength_plot.RData")
