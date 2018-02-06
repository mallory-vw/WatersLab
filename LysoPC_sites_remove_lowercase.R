#####load workspace#####
load("J:/III/Waters/Group Members/Mallory/R/LysoPC_sites_remove_lowercase.RData")

#####setwd#####
setwd("J:/III/Waters/Group Members/Mallory/LysoPC_project/FIRE_results/Motifs/MotifSitesFiles/")


#####load data#####
filenames <- dir()
filenames <- filenames[c(1,10:17,2:9)]
motif_data <- lapply(filenames, read.csv, stringsAsFactors = FALSE, header=F)
names(motif_data) <- filenames

#####remove lowercase letters#####
motif_data_onlySite <- lapply(motif_data, function(j){
  for (i in seq_along(j$V1)){
    j[i,1] <- gsub('[a-z]', '', j[i,1])
  }
  j
})

filenames_onlySite <- gsub('\\.','_onlySite.', filenames)
names(motif_data_onlySite) <- filenames_onlySite

#just pull them out so viewing is easier, don't need to do this part
list2env(motif_data, envir=.GlobalEnv)
list2env(motif_data_onlySite, envir=.GlobalEnv)


#####write csvs#####
#can lapply the list of only site, will automatically make the text files
#wanted it saved with just a LF end-of-line (unix format) - argument eol in write.table
#needed to open the file in binary first?  this is the line that assigns f, just use exactly what you would use for the filename in the 
#write.table function, but also add open = "wb" - i think this opens the connection to the file as binary, which will allow me 
#to set the unix end-of-line without windows interfering
#then use the normal write.table, and set file = f, which was already defined in the function
# note: f will change for each item in motif_data_onlySite, so will be reset fo each motif the same way it would in a forloop
#set eol = "\n", which should just set it to LF (linefeed), the unix default ending
#then close f, the connection to the file
#http://r.789695.n4.nabble.com/Write-table-eol-argument-td3159733.html

lapply(1:length(motif_data_onlySite), function(q){
  f <- file(paste(names(motif_data_onlySite[q]),sep=""),open="wb")
  write.table(motif_data_onlySite[[q]],
            file = f,
            col.names=FALSE,
            row.names = FALSE,
            quote = FALSE,
            eol = "\n")
  close(f)
})




#####save workspace#####
save.image("J:/III/Waters/Group Members/Mallory/R/LysoPC_sites_remove_lowercase.RData")
