#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/Brett_GeneListIntersect.RData")

#####setwd######
setwd("J:/III/Waters/Group Members/Mallory/Helping/Brett_GeneListIntersect/")

# install.packages('VennDiagram')
library(VennDiagram)

#####load data#####
Control <- read.csv("AllControl.csv", header = F)[,1]
Exp1 <- read.csv("Exp1.csv", header = F)[,1]
Exp2 <- read.csv("Exp2.csv", header = F)[,1]
Exp3 <- read.csv("Exp3.csv", header = F)[,1]


####ControlIncluded####
####In Exp1, not in Exp2
Exp1_NotExp2 <- setdiff(Exp1, Exp2)
  write.csv(Exp1_NotExp2,"Exp1_NotExp2.csv")
####In Exp1, not in Exp3
Exp1_NotExp3 <- setdiff(Exp1, Exp3)
  write.csv(Exp1_NotExp3,"Exp1_NotExp3.csv")

####In Exp2, not in Exp1
Exp2_NotExp1 <- setdiff(Exp2, Exp1)
  write.csv(Exp2_NotExp1,"Exp2_NotExp1.csv")
####In Exp2, not in Exp3
Exp2_NotExp3 <- setdiff(Exp2, Exp3)
  write.csv(Exp2_NotExp3,"Exp2_NotExp3.csv")

####In Exp3, not in Exp1
Exp3_NotExp1 <- setdiff(Exp3, Exp1)
  write.csv(Exp3_NotExp1,"Exp3_NotExp1.csv")
####In Exp3, not in Exp2
Exp3_NotExp2 <- setdiff(Exp3, Exp2)
  write.csv(Exp3_NotExp2,"Exp3_NotExp2.csv")
  
####InExp1, Not in Exp2+Exp3  
Exp1_NotExp2_Exp3 <- setdiff(Exp1, union(Exp2,Exp3))
  write.csv(Exp1_NotExp2_Exp3,"Exp1_NotExp2_Exp3.csv")
####InExp2, Not in Exp1+Exp3  
Exp2_NotExp1_Exp3 <- setdiff(Exp2, union(Exp1,Exp3))
  write.csv(Exp2_NotExp1_Exp3,"Exp2_NotExp1_Exp3.csv")
####InExp3, Not in Exp2+Exp1
Exp3_NotExp2_Exp1 <- setdiff(Exp3, union(Exp2,Exp1))
  write.csv(Exp3_NotExp2_Exp1,"Exp3_NotExp2_Exp1.csv")
  
  
  
##Shared
Exp1_Exp2 <- intersect(Exp1, Exp2)
  write.csv(Exp1_Exp2,"Exp1_Exp2.csv")
Exp1_Exp3 <- intersect(Exp1, Exp3)
  write.csv(Exp1_Exp3,"Exp1_Exp3.csv")
Exp2_Exp3 <- intersect(Exp2, Exp3)  
  write.csv(Exp2_Exp3,"Exp2_Exp3.csv")
Exp1_Exp2_Exp3 <- Reduce(intersect, list(Exp1,Exp2,Exp3))
  write.csv(Exp1_Exp2_Exp3,"Exp1_Exp2_Exp3.csv")

####In Control, not in Exp1
Control_NotExp1 <- setdiff(Control, Exp1)
  write.csv(Control_NotExp1, "Control_NotExp1.csv")
####In Control, not in Exp2
Control_NotExp2 <- setdiff(Control, Exp2)
  write.csv(Control_NotExp2, "Control_NotExp2.csv")
####In Control, not in Exp3
Control_NotExp3 <- setdiff(Control, Exp3)
  write.csv(Control_NotExp3, "Control_NotExp3.csv")

####In Control, In Exp1
Control_Exp1 <- intersect(Control,Exp1)
  write.csv(Control_Exp1,"Control_Exp1.csv")
####In Control, In Exp2
Control_Exp2 <- intersect(Control,Exp2)
  write.csv(Control_Exp2,"Control_Exp2.csv")
####In Control, In Exp3
Control_Exp3 <- intersect(Control,Exp3)
  write.csv(Control_Exp3,"Control_Exp3.csv")
  
#####ControlExcluded#####
####In Exp1, not in Control
Exp1_NotControl <- setdiff(Exp1, Control)
  write.csv(Exp1_NotControl, "Exp1_NotControl.csv")
####In Exp2, not in Control
Exp2_NotControl <- setdiff(Exp2, Control)
  write.csv(Exp2_NotControl, "Exp2_NotControl.csv")
####In Exp3, not in Contol
Exp3_NotControl <- setdiff(Exp3, Control)
  write.csv(Exp3_NotControl, "Exp3_NotControl.csv")

####In Exp1_NotControl, not in Exp2_NotControl
Exp1_Unique_NotExp2_unique <- setdiff(Exp1_NotControl, Exp2_NotControl)
  write.csv(Exp1_Unique_NotExp2_unique, "Exp1_Unique_NotExp2_unique.csv")
####In Exp1_NotControl, not in Exp3_NotControl
Exp1_Unique_NotExp3_unique <- setdiff(Exp1_NotControl, Exp3_NotControl)
  write.csv(Exp1_Unique_NotExp3_unique, "Exp1_Unique_NotExp3_unique.csv")

####In Exp2_NotControl, not in Exp1_NotControl
Exp2_Unique_NotExp1_unique <- setdiff(Exp2_NotControl, Exp1_NotControl)
  write.csv(Exp2_Unique_NotExp1_unique, "Exp2_Unique_NotExp1_unique.csv")
####In Exp2_NotControl, not in Exp3_NotControl
Exp2_Unique_NotExp3_unique <- setdiff(Exp2_NotControl, Exp3_NotControl)
  write.csv(Exp2_Unique_NotExp3_unique, "Exp2_Unique_NotExp3_unique.csv")

####In Exp3_NotControl, not in Exp1_NotControl
Exp3_Unique_NotExp1_unique <- setdiff(Exp3_NotControl, Exp1_NotControl)
  write.csv(Exp3_Unique_NotExp1_unique, "Exp3_Unique_NotExp1_unique.csv")
####In Exp3_NotControl, not in Exp2_NotControl
Exp3_Unique_NotExp2_unique <- setdiff(Exp3_NotControl, Exp2_NotControl)
  write.csv(Exp3_Unique_NotExp2_unique, "Exp3_Unique_NotExp2_unique.csv")

####Shared
Exp1_NotControl_Exp2_NotControl <- intersect(Exp1_NotControl, Exp2_NotControl)
  write.csv(Exp1_NotControl_Exp2_NotControl, "Exp1_NotControl_Exp2_NotControl.csv")
Exp1_NotControl_Exp3_NotControl <- intersect(Exp1_NotControl, Exp3_NotControl)
  write.csv(Exp1_NotControl_Exp3_NotControl, "Exp1_NotControl_Exp3_NotControl.csv")
Exp2_NotControl_Exp3_NotControl <- intersect(Exp2_NotControl, Exp3_NotControl)
  write.csv(Exp2_NotControl_Exp3_NotControl, "Exp2_NotControl_Exp3_NotControl.csv")


#####CompletelyUnique#####
Control_completely_Unique <- Reduce(setdiff, list(Control,Exp1,Exp2,Exp3))
  write.csv(Control_completely_Unique, "Control_completely_Unique.csv")
Exp1_completely_Unique <- Reduce(setdiff, list(Exp1,Control,Exp2,Exp3))
  write.csv(Exp1_completely_Unique, "Exp1_completely_Unique.csv")
Exp2_completely_Unique <- Reduce(setdiff, list(Exp2,Control,Exp1,Exp3))
  write.csv(Exp2_completely_Unique, "Exp2_completely_Unique.csv")
Exp3_completely_Unique <- Reduce(setdiff, list(Exp3,Control,Exp2,Exp1))
  write.csv(Exp3_completely_Unique, "Exp3_completely_Unique.csv")

####venn diagram####
#area1 - 86 - control
#area2 - 21 - exp1
#area3 - 80 - exp2
#area4 - 65 - exp3
#n12
length(intersect(Control,Exp1))
  #13
#n13
length(intersect(Control, Exp2))
  #39
#n14
length(intersect(Control, Exp3))
  #46
#n23
length(intersect(Exp1,Exp2))
  #9
#n24
length(intersect(Exp1,Exp3))
  #9
#n34
length(intersect(Exp2,Exp3))
  #39
#n123
length(intersect(intersect(Control, Exp1), Exp2))
  #7
#n124 
length(intersect(intersect(Control, Exp1), Exp3))
  #9
#n134
length(intersect(intersect(Control, Exp2), Exp3))
  #30
#n234
length(intersect(intersect(Exp1, Exp2), Exp3))
  #6
#n1234
length(intersect(intersect(intersect(Control,Exp1),Exp2), Exp3))
  #6

png("AllGenes_venndiagram.png", 490, 490)
draw.quad.venn(area1 = 86, area2 = 21, area3 = 80, area4 = 65, 
               n12 = 13, n13 = 39, n14 = 46, 
               n23 = 9, n24 = 9, 
               n34 = 39, 
               n123= 7, n124 = 9, n134 = 30, n234 = 6,
               n1234 = 6,
                 category=(c("Control","Exp1","Exp2","Exp3")),
                 fill = c("red","blue","yellow","green"))
dev.off()


####just the NotControl lists
#area1 - Exp1_NotContol - 8
#area2 - Exp2_NotControl - 41
# area3 - Exp3_NotControl - 19

#n12
length(intersect(Exp1_NotControl,Exp2_NotControl))
#2

#n13
length(intersect(Exp1_NotControl,Exp3_NotControl))
#0

#n23
length(intersect(Exp2_NotControl,Exp3_NotControl))
#9

#n123
length(intersect(intersect(Exp1_NotControl, Exp2_NotControl), Exp3_NotControl))
#0


png("NotControl_venndiagram.png", 490, 490)
draw.triple.venn(area1 = 8, area2 = 41, area3 = 19, 
                                             n12 = 2, n13 = 0, n23 = 9, 
                                             n123= 0, 
                                             category=(c("Exp1","Exp2","Exp3")),
                                             fill = c("blue","yellow","green"))
dev.off()


#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/Brett_GeneListIntersect.RData")
