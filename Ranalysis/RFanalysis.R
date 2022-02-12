#Trees:

#55 Divergence Times:
#	12 Exons
#		8 recombination rates
#			5 lengths (1024)
#				7 mutation scales u_.15.ali or u_1.5.ali
#					*100 seqgen reps
#						*100 reps

#Take bootstrap trees generated for some recombination value at 1024 length
#Find the msout trees they were generated from
#perform RF analysis on msout trees <-> seqgen_Raxml trees. For each exon

#10 seqgen trees with 100 raxml bootstrap trees compared to msout trees

#cd T_80000_40000/r0/1/1024

#head ../EXONS_100trees.tre 
#head u_15.ali
#head RAxML_bootstrap.u_15_1


#RF distance of original 10 seqgen alignments vs each alignment's 100 raxml bootstrap replicates
filepath <- "/pool/Recombination-Michael/6_29_18/"
T1_T2 <- "T_80000_10000"
pdf(paste(T1_T2, "_all.pdf", sep=""))
r_rate <- "r0"
RFdist <- array(0, dim=c(1000,8))
ms_trees = paste(filepath,T1_T2,"/1/",r_rate,"/EXONS_100trees.tre", sep="")
ms_data <- read.table(ms_trees, sep=" ", stringsAsFactors=FALSE)
scale_count = 0

for(scale in c("1.5", ".15", ".015", ".0015", ".00015", ".000015", ".0000015")){
	scale_count=scale_count + 1
	count=1
	for(i in seq(1,10)){
		t1 <- read.tree(text=ms_data[i,1])		#Get the simulated tree from 100Exons

		raxml_trees = paste(filepath,"reduced_parameter_sets/",T1_T2,"/1/", r_rate, "/1024/u_",scale,"/RAxML_bootstrap.u_", scale, "_", i, sep="")
		raxml_data <- read.table(raxml_trees, sep=" ", stringsAsFactors=FALSE) #get the corresponding raxml estimation based on the right seqgen rep

		RF_total = 0

                for(j in seq(1,100)){
                        t2 <- read.tree(text=raxml_data[j,1])
                        RF_total = RF_total + dist.topo(t1,t2, method="PH85")
                }
                RFdist[count,scale_count] <- RF_total/100
                count=count+1
	}
}

plot(0,0, xlim=c(0,10), ylim=c(0,45))
points(rep(1,1000)+runif(1000,-0.4,0.4), RFdist[,1], cex=0.2)
points(rep(2,1000)+runif(1000,-0.4,0.4), RFdist[,2], cex=0.2)
points(rep(3,1000)+runif(1000,-0.4,0.4), RFdist[,3], cex=0.2)
points(rep(4,1000)+runif(1000,-0.4,0.4), RFdist[,4], cex=0.2)
points(rep(5,1000)+runif(1000,-0.4,0.4), RFdist[,5], cex=0.2)
points(rep(6,1000)+runif(1000,-0.4,0.4), RFdist[,6], cex=0.2)
points(rep(7,1000)+runif(1000,-0.4,0.4), RFdist[,7], cex=0.2)
dev.off()


####For average RF dist across 100 bootstrap replicates


#plot(0,0, xlim=c(0,10), ylim=c(0,45), ylab="RF distance")

library(ape)
filepath <- "/pool/Recombination-Michael/6_29_18/"
color_space <- c("red","blue","green")

for(T1_T2 in c("T_80000_10000","T_1280000_2500","T_2560000_10000","T_2560000_1280000","T_2560000_20000","T_320000_10000","T_320000_20000","T_5000_2500","T_640000_320000","T_80000_10000","T_80000_5000")){
#for(T1_T2 in c("T_80000_10000")){
	#pdf(paste(T1_T2, "_all_weighted.pdf", sep=""))
	pdf(paste(T1_T2, "_all_unrooted.pdf", sep=""))
	plot(0,0,col="white", xaxt="n",xlim=c(0,8), ylim=c(45,0), ylab="RF distance", xlab="Substitution Rate")
	#T1_T2 <- "T_80000_10000"
	color_iter = 0
	RFmeans <- array(0, dim=c(3,7))
	lengths = 1
	for(len in c("64","256","1024")){
		color_iter = color_iter + 1;
		color <- color_space[color_iter]
		r_rate <- "r0"
		RFdist <- array(0, dim=c(100,7))
		dimnames(RFdist)[[2]] <- c("1.5", ".15", ".015", ".0015", ".00015", ".000015", ".0000015")
		ms_trees = paste(filepath,T1_T2,"/1/",r_rate,"/EXONS_100trees.tre", sep="")
		ms_data <- read.table(ms_trees, sep=" ", stringsAsFactors=FALSE)
		scale_count = 0
		for(scale in c("1.5", ".15", ".015", ".0015", ".00015", ".000015", ".0000015")){
			scale_count=scale_count + 1
			count=1
			for(i in seq(1,100)){
				t1 <- read.tree(text=ms_data[i,1])		#Get the simulated tree from 100Exons
				unrooted_t1 <- unroot(t1)
				raxml_trees = paste(filepath,"reduced_parameter_sets/",T1_T2,"/1/",r_rate,"/",len,"/u_",scale,"/RAxML_bootstrap.u_", scale, "_", i, sep="")
				raxml_data <- read.table(raxml_trees, sep=" ", stringsAsFactors=FALSE) #get the corresponding raxml estimation based on the right seqgen rep
				RF_total = 0
				for(j in seq(1,100)){
					t2 <- read.tree(text=raxml_data[j,1])
					unrooted_t2 <- unroot(t2)
					#RF_total = RF_total + dist.topo(t1,t2, method="score")
					#RF_total = RF_total + dist.topo(t1,t2, method="PH85")
					RF_total = RF_total + dist.topo(unrooted_t1,unrooted_t2, method="PH85")
				}
				RFdist[count,scale_count] <- (RF_total/100)
				count=count+1
			}
		}
		RFmeans[lengths,] <- colMeans(RFdist)
		#plot(0,0, xlim=c(0,100), ylim=c(0,45))
		points(rep(1,100)+runif(100,-0.4,0.4), RFdist[,1], cex=0.5, col=color)
		points(rep(2,100)+runif(100,-0.4,0.4), RFdist[,2], cex=0.5, col=color)
		points(rep(3,100)+runif(100,-0.4,0.4), RFdist[,3], cex=0.5, col=color)
		points(rep(4,100)+runif(100,-0.4,0.4), RFdist[,4], cex=0.5, col=color)
		points(rep(5,100)+runif(100,-0.4,0.4), RFdist[,5], cex=0.5, col=color)
		points(rep(6,100)+runif(100,-0.4,0.4), RFdist[,6], cex=0.5, col=color)
		points(rep(7,100)+runif(100,-0.4,0.4), RFdist[,7], cex=0.5, col=color)
		lengths = lengths+1
	}
	lengths = 1;
	for(color in color_space){
		lines(RFmeans[lengths,], col=color, lwd=3)
	        points(RFmeans[lengths,], pch=24, col="black", bg=color, cex=1.5)
		lengths = lengths+1
	}
	title(T1_T2)
	legend(6,20, c("64","256","1024"), pch=2, col=color_space, title="Locus Length")
	par(cex.axis=.8)
	axis(1, at=1:7, labels=c("1.5", ".15", ".015", ".0015", ".00015", ".000015", ".0000015"))
	dev.off()
	#dev.print(pdf, paste(T1_T2, "_weighted.pdf", sep=""))
}
