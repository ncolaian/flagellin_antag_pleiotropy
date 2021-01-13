#go through and perform mit analyisis on results that were trimmed by trimal

#MIT analysis for the 

library(MASS)
rcw_conservation <- function(x) {
  new <- x
  for ( i in 1:nrow(x) ) {
    pos1 <- x[i,1]
    pos2 <- x[i,2]
    
    
    sum1 <- sum(x$Entropy[x$Pos1 == pos1 | x$Pos2 == pos1])
    sum2 <- sum(x$Entropy[x$Pos1 == pos2 | x$Pos2 == pos2])
    #Perform the calculation
    new[i,3] <- x[i,3]/
      (( (sum1+sum2) - (2*x[i,3]) )/
         ( length(x$Entropy[x$Pos1 == pos2 | x$Pos2 == pos2]) + length(x$Entropy[x$Pos1 == pos1 | x$Pos2 == pos1]) - 2 ))
  }
  return(new)
}

z_score <- function(vect_numbers) {
  out <- c()
  for (i in vect_numbers) {
    out <- append(out, (i - mean(vect_numbers))/sd(vect_numbers) )
  }
  return(out)
}

run_analysis <- function(mit, msa = "NA") {
  
  if ( msa != "NA" ) {
    msa$pos <- 1:length(msa[,1])
    mit <- mit[msa$V1[mit$Pos1] > 7 & msa$V1[mit$Pos2] > 7,]
  }
  
  #row column weighting
  mit <- rcw_conservation(mit)
  mit <- na.omit(mit)
  
  #get z-score
  mit$zscore <- z_score(mit$Entropy)
  
  #assign p-vals
  mit$pvalue <- pnorm(mit$zscore, lower.tail = F)
  
  return(mit)
}

###########
## FILES ##
###########
# mit file from the karasov genomes
mit_file <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/karasov_fliC/mit.txt")

##################################################################
## ONE FUNCTION TO perform RCW, Z-SCORE, AND PVALUE CALCULATION ##
##################################################################
mit_file <- run_analysis(mit_file)

#Read in the actual location of each position on the fliC protein
col_position <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/karasov_fliC/correct_position.txt", header = T)

col_position$Pseudomonas_aeruginosa_PASGNDM699_7104 <- col_position$Pseudomonas_syringae_pv_tomato_DC3000_111
col_position$Pseudomonas_syringae_pv_tomato_DC3000_111 <- row.names(col_position)
mit_file$plot_syringe_1 <- col_position$Pseudomonas_syringae_pv_tomato_DC3000_111[mit_file$Pos1]
mit_file$plot_syringe_2 <- col_position$Pseudomonas_syringae_pv_tomato_DC3000_111[mit_file$Pos2]

mit_file$plot_aerug_1 <- col_position$Pseudomonas_aeruginosa_PASGNDM699_7104[mit_file$Pos1]
mit_file$plot_aerug_2 <- col_position$Pseudomonas_aeruginosa_PASGNDM699_7104[mit_file$Pos2]


mean(mit_file$Entropy)
sd(mit_file$Entropy)
quantile(mit_file$Entropy, .5)
hist(mit_file$Entropy)
hist(mit_file$zscore) 

#get p-values from multiple testing with bonferroni
top <- mit_file[mit_file$Entropy > quantile(mit_file$Entropy,.995),]
#write.csv(top,"/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj//882_data/updated_work_10_23_18/new_clusters/mit/clade1_top_hits.csv", row.names = F)
write.csv(top,"/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/karasov_fliC/top_mit_positions.csv", row.names = F)


#display the results
#this is really just getting the most significant z_scores
#top_combined <- c(top$plot_Pos1, top$plot_Pos2)
top_combined <- c(top$Pos1, top$Pos2)
sort(unique(top_combined))
top_combined <- na.omit(top_combined)


#graph
bonf_trial <- hist(top_combined, breaks = (max(top_combined)-min(top_combined)), ylim = c(0,20), main = "Top MIT results by position",
                   xlab = "Position", ylab = "# of Results in the top 0.5%", xaxt = "n")+
  axis(1, c(1,30,43,51, 100, 150,200,250 ), labels = c(1,30,43,51, 100, 150,200,250 ), tick = T,padj= -1)

library(ggplot2)

#try a manhatten plot
#library(qqman)
#manhattan(mit_file, chr = "plot_Pos1", bp = "plot_Pos2", p= "pvalue", snp = "Entropy")

#normal
mit_file$plot_Pos1 <- factor(as.character(mit_file$Pos1), levels = as.character(1:max(mit_file$Pos1)) )
mit_file$plot_Pos2 <- factor(as.character(mit_file$Pos2), levels = as.character(1:max(mit_file$Pos2)) )

mit_file$plot_Entropy <- mit_file$Entropy
mit_file$plot_Entropy[mit_file$plot_Entropy < 0] <- 0
mit_file$plot_Entropy[mit_file$plot_Entropy < quantile(mit_file$plot_Entropy,.995)] <- 0

flg22_pos <- c(30,43,51)
ggplot(mit_file, aes(x=plot_Pos2, y=plot_Pos1))+
  geom_tile(aes(fill=plot_Entropy), colour="white",size=.0001)+
  scale_fill_gradient(high = "firebrick", low = "white")+
  scale_y_discrete(expand=c(0,0),position = "right", breaks=flg22_pos)+
  scale_x_discrete(expand=c(0,0), breaks=flg22_pos)+
  #scale_y_discrete(expand=c(0,0),position = "right", breaks=c(seq(col_position$X1[1],max(col_position$X1), by = 100), flg22_pos, max(col_position$X1)-1))+
  #scale_x_discrete(expand=c(0,0), breaks=c(seq(col_position$X1[2],max(col_position$X1), by = 100), flg22_pos, max(col_position$X1)))+
  coord_fixed()+
  xlab("Position 1")+
  ylab("Position 2")+
  labs(fill="MIT Entropy")+
  theme_dark()

#I want to show that the flg22 region is overrepresented in MIT positions indicating this region is 
#differentiating with other positions. It is no fluke that there are so many positions in the flg22 region in the top MIT
#I will bootstrap.

set.seed(1994)

positions <- unique(c(mit_file$Pos1, mit_file$Pos2))
flg22_positions <- seq(30, 51, by=1)

number_of_pulls <- nrow(top)
#there are about 200 top MIT results

number_of_flg22_seqs <- sum(top$Pos1 %in% flg22_positions | top$Pos2 %in% flg22_positions)

#Go ahead and see how likely it is to have more in this region than the results above
#10000 draws
how_many_contained_flg22 <- c()
for (i in 1:10000) {
  count <- 0
  for ( j in 1:number_of_pulls) {
    if (sum( sample(positions, 2, replace = F) %in% flg22_positions ) > 0 ) {
      count <- count+1
    }
  }
  how_many_contained_flg22 <- append(how_many_contained_flg22, count)
}

#How many random draws had equal or more positions in the flg22 region?
sum(how_many_contained_flg22 >= number_of_flg22_seqs)/10000
max(how_many_contained_flg22)
hist(how_many_contained_flg22)

#the flg22 region is enriched in MIT regions

#Now I need to look at positive or negative selection
#approximately half of the high MIT positions occur in the flg22 region
