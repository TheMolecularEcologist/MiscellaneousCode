#############################################################################
### Script for viewing and manipulating tabular output from BLAST searches
### Use the (-m 8) otpion to get a tabular output from blastall
### Not very fast at the moment, may need som optimisation
#############################################################################
### Robert Ekblom
### 2/3-2009
### Thanks to Owen Petchey for help with the loop
#############################################################################

rm(list=ls())   #Emties the memory
setwd("C:\\blast\\")  #Setting the working directory

### Open the file ###
fil <- file.choose()
out <- read.table(fil, header=FALSE)

### Add informative column headings ###
names(out)<- c("Query_id",
		"Subject_id",
		"%_identity",
		"length",
		"mismatches",
		"gaps",
		"q_start",
		"q_end",
		"s_start",
		"s_end",
		"e_value",
		"score")

### Add direction column ###
out$strand <- sign(out$s_end-out$s_start)
# attach(out) #make variables acssible by name TURNED OFF!!!

### Sorting function, Not needed fo the rest of the script! ###
#out.sort1 <- with(out,out[order(Query_id, Subject_id, -strand, s_start),  ])
head(out.sort1)
str(out.sort1)

out <- with(out,out[order(Query_id, e_value, -score, -length),  ])
head(out)

### Take only the best match for each query

all.genes <- unique(out$Query_id)
for(i in 1:length(all.genes)){
		goi <- out$Query_id == all.genes[i]
		answ <- out[goi,][which(out$e_value[goi] == min(out$e_value[goi]))[1], ]
		if(i == 1)
			out.new <- answ
		if(i > 1)
			out.new <- rbind(out.new, answ)
		}

### Take only rows with an e-value less than "maxe"
maxe <- 1e-5  #Change here to optimise
out.new <- out.new[out.new$e_value < maxe,]

### Get some summary info
head(out.new)
str(out.new)
paste("Name of BLASTA output file: ", fil)
paste("Number of unique queries with hits: ", length(all.genes))
paste("Number of unique hits: ", length(unique(out$Subject_id)))
paste("Cut-off e-value value: ", maxe)
paste("Number of unique queries below cut-off e: ", length(out.new$Query_id))

### Export tables 
write.table(out.new, "output\\yourfile.rout", quote=FALSE, row.names=FALSE)
write.table(out.new$Query_id, "output\\yourfile.queryid.rout", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(out.new$Subject_id, "output\\yourfile.subjectid.rout", quote=FALSE, row.names=FALSE, col.names=FALSE)
