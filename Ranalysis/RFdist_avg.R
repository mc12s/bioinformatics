#!/usr/bin/Rscript
#Rscript RFdist.R 5.0E-6 12
library(ape)
library(stringr)
library(fields)
library(colorspace)

args <- commandArgs(trailingOnly = TRUE)
#rrate <- args[1]
#reps <- as.numeric(args[2])

pal <- choose_palette()

recom_array <- matrix(NA,11,11)
recom_array <- array(c(recom_array,recom_array),dim=c(11,11,9))

for(r_iter in seq(1:6)){
#for(r_iter in seq(1:1)){

#rrates <- c("0.5","0.005","5.0E-4","5.0E-5","5.0E-6","5.0E-7","0.0")
rrates <- c(".02",".2","2","20","200","0")
#rrates <- c("0.5")
rrate <- rrates[r_iter]
reps <- 10000
print(reps)
#fileprefix <- "T_5000_2500/tree_10000_97500_2_0.0__"


#heatmapdata <- array(NA, dim=c(32,8,8))

files = c("T_5000_2500/",
					"T_10000_2500/",
					"T_20000_2500/",
					"T_40000_2500/",
					"T_80000_2500/",
					"T_160000_2500/",
					"T_320000_2500/",
					"T_640000_2500/",
					"T_1280000_2500/",
					"T_2560000_2500/",
					"T_10000_5000/",
					"T_20000_5000/",
					"T_40000_5000/",
					"T_80000_5000/",
					"T_160000_5000/",
					"T_320000_5000/",
					"T_640000_5000/",
					"T_1280000_5000/",
					"T_2560000_5000/",
					"T_20000_10000/",
					"T_40000_10000/",
					"T_80000_10000/",
					"T_160000_10000/",
					"T_320000_10000/",
					"T_640000_10000/",
					"T_1280000_10000/",
					"T_2560000_10000/",
					"T_40000_20000/",
					"T_80000_20000/",
					"T_160000_20000/",
					"T_320000_20000/",
					"T_640000_20000/",
					"T_1280000_20000/",
					"T_2560000_20000/",
					"T_80000_40000/",
					"T_160000_40000/",
					"T_320000_40000/",
					"T_640000_40000/",
					"T_1280000_40000/",
					"T_2560000_40000/",
					"T_160000_80000/",
					"T_320000_80000/",
					"T_640000_80000/",
					"T_1280000_80000/",
					"T_2560000_80000/",
					"T_320000_160000/",
					"T_640000_160000/",
					"T_1280000_160000/",
					"T_2560000_160000/",
					"T_640000_320000/",
					"T_1280000_320000/",
					"T_2560000_320000/",
					"T_1280000_640000/",
					"T_2560000_640000/",
					"T_2560000_1280000/"
					)

heatmapRF.data <- data.frame(matrix(NA, 32, length(files)))
heatmapScores.data <- data.frame(matrix(NA, 32, length(files)))

#files <- paste(files,sep="")

colname <- files
colname <- str_split_fixed(colname,"/",2)[,1]
colname <- str_split_fixed(colname,"T_",2)[,2]
heatmapRF.data <- setNames(heatmapRF.data, colname)
heatmapScores.data <- setNames(heatmapScores.data, colname)
RFcolsum.data <- data.frame(matrix(NA, reps, 28))
scorescolsum.data <- data.frame(matrix(NA, reps, 28))

mapcol <- 1


recom_array <- matrix(NA,11,11)



#lengths = c("12")
#length_ind=1
#for (len in lengths){
RFavg <- array(c(NA), dim = c(55,6))
len <- 12
recom_rates = c(".02",".2","2","20","200","20000")
#rate_ind=1
	for (rate in recom_rates){
		#RFind = 1

		rate <- "2"
		len <- 12
		for(fileprefix in files){

			treefile <- paste(fileprefix,"10000reps_",rate,"_",len,".msoutextracted1_2.txt",sep="")
			treedata <- read.table(treefile, sep=" ", stringsAsFactors=FALSE)
			RFdist <- 0
			for(rep in seq(1,19999,by=2)){
				if((rep%%100)+1==0){
					print(rep)
				}
				spl <- regexpr("\\(",treedata[rep,])
				t1 <- substring(treedata[rep,],c(spl,nchar(treedata[rep,])))[1]
				t1 <- read.tree(text=t1)
				spl <- regexpr("\\(",treedata[rep+1,])
				t2 <- substring(treedata[rep+1,],c(spl,nchar(treedata[rep+1,])))[1]
				t2 <- read.tree(text=t2)
				RFdist <- RFdist + dist.topo(t1,t2, method="PH85")
			}
			print(treefile)
			print(RFdist/10000)
			#RFind=RFind+1

		}
		rate_ind=rate_ind+1
	}
#	length_ind = length_ind+1
#}

	repcol <- 1
	template_matrix <- matrix(0,32,reps+1)
	allRFdist.data <- data.frame(template_matrix)
	allscores.data <- data.frame(template_matrix)
	allRFdist.data[1] <- seq(1,64,by=2)
	allscores.data[1] <- seq(1,64,by=2)
	for(rep in seq(1,reps)){
		repcol <- repcol+1
		#if(r_iter==1){
		treefile <- paste(fileprefix,"10000reps_.2_12.msoutEXONS.out",sep="")
		#}else{
		#	treefile <- paste(fileprefix,"__",(rep+90),".tre",sep="")
		#}
		
		alldist.data <- data.frame(matrix(0,32,2))
	
		trytree <- function(treefile){
			out <- tryCatch(
				{
					print(treefile)
					if(file.exists(treefile)){ #limit the broken connection
						treedata <- read.table(treefile, sep=" ", stringsAsFactors=FALSE)
					}else{
						alldist.data[1] <- data.frame(rep(NA, 32))
						alldist.data[2] <- data.frame(rep(NA, 32))
						return(alldist.data)
					}
					if(treedata[1,1]=="n8_c0_a0_Did_not_Coalesce" || treedata[1,1]=="n8_cx_a0_Did_not_Coalesce" || treedata[1,1]=="n8_cX_a0_Did_not_Coalesce" ){			
						alldist.data[1] <- data.frame(rep(NA, 32))
						alldist.data[2] <- data.frame(rep(NA, 32))
						return(alldist.data)
					}
					RFdist <- vector(mode="integer", length=32)
					branchscores <- vector(mode="double", length=32)

					counter=1
					for(i in seq(1,64,by=2)){
						t1 <- read.tree(text=treedata[i,2])
						#plot(t1)
						#edgelabels(t1$edge.length, bg="black", col="white", font=2)
						t2 <- read.tree(text=treedata[i+1,2])
						RFdist[counter] <- dist.topo(t1,t2, method="PH85")
						##branchscores[counter] <- dist.topo(t1,t2, method="score")
						counter = counter+1
					}
					alldist.data[1] <- data.frame(RFdist)
					alldist.data[2] <- data.frame(branchscores)
					return(alldist.data)
				},
				error=function(cond) {
					#message(paste("File error: ",treefile))
					message(cond)
					alldist.data[1] <- data.frame(rep(NA, 32))
					alldist.data[2] <- data.frame(rep(NA, 32))
					return(alldist.data)
				},
				warning=function(cond){
					message(cond)
					alldist.data[1] <- data.frame(rep(NA, 32))
					alldist.data[2] <- data.frame(rep(NA, 32))
					return(alldist.data)
				}
			)
			return(out)
		}

		alldist.data <- trytree(treefile)
		#print(alldist.data)
		allRFdist.data[repcol] <- alldist.data[1]
		#print(allRFdist.data)
		##allscores.data[repcol] <- alldist.data[2]
	
	}


	RFoutfile <- paste(unlist(strsplit(fileprefix,"/"))[1],"/RFtopo_",unlist(strsplit(fileprefix,"/"))[2],".score", sep="")
	print(RFoutfile)
	write.table(allRFdist.data, sep="\t", file=RFoutfile, row.names=FALSE)
	print(RFoutfile)
	heatmapRF.data[mapcol] <- rowMeans(allRFdist.data[,-1], na.rm=TRUE)

	RFcolsum.data[mapcol] <- colSums(allRFdist.data[,-1], na.rm=TRUE)
	

	##scoresOutfile <- paste(unlist(strsplit(fileprefix,"/"))[1],"/RFbranch_",unlist(strsplit(fileprefix,"/"))[2],".score", sep="")
	##write.table(allscores.data, sep="\t", file=scoresOutfile, row.names=FALSE)
	##print(scoresOutfile)
	##heatmapScores.data[mapcol] <- rowMeans(allscores.data[,-1], na.rm=TRUE)
	##colSums(allscores.data[,-1], na.rm=TRUE, dims=1)

	##scorescolsum.data[mapcol] <- colSums(allscores.data[,-1], na.rm=TRUE)

	
	mapcol <- mapcol+1
}

print(rowSums(RFcolsum.data))
print(rowSums(scorescolsum.data))

#print(heatmapRF.data)
#print(heatmapScores.data)

#allheatmaps <- matrix(0,11,11)
#iter <-1
#for(i in seq(3,32,by=4)){
#	cur_heatmap<- matrix(0,11,11)
	#cur_heatmap[lower.tri(cur_heatmap,diag=FALSE)] <- as.matrix(heatmapRF.data)[i,]
	#cur_heatmap <- t(cur_heatmap)
#	cur_heatmap[lower.tri(cur_heatmap,diag=FALSE)] <- as.matrix(heatmapScores.data)[i,]
#	allheatmaps <- array(c(allheatmaps,cur_heatmap), dim=c(11,11,iter+1)) 
	#image(d)
#	iter <- iter+1
#}

#axis_names <- c("2500","5000","10000","20000","40000","80000","160000","320000","640000","1280000","2560000")
#dimnames(allheatmaps) <- list(axis_names, axis_names)

#par(mfrow=c(8,5))
#par(mar=c(5,5,5,5))
#for(i in seq(2,9)){
#	image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
#	axis( 2, at=seq(0,1,length.out=nrow( allheatmaps[,,i] ) ), labels= rownames( allheatmaps[,,i] ), las= 2 )
#	axis( 1, at=seq(0,1,length.out=ncol( allheatmaps[,,i] ) ), labels= colnames( allheatmaps[,,i] ), las= 2 )
#}

#par(mar=c(0,0,0,0))
#par(mfrow=c(8,1))
#for(i in seq(2,9)){
#	image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
#	box()
#}
#dev.print(pdf,paste(rrate,'r6heat.pdf', sep=""))

allheatmaps <- matrix(NA,11,11)
iter <-1
for(i in seq(3,32,by=4)){
	cur_heatmap<- matrix(NA,11,11)
	cur_heatmap[lower.tri(cur_heatmap,diag=FALSE)] <- as.matrix(heatmapRF.data)[i,]
	#cur_heatmap <- t(cur_heatmap)
	#cur_heatmap[lower.tri(cur_heatmap,diag=FALSE)] <- as.matrix(heatmapScores.data)[i,]
	allheatmaps <- array(c(allheatmaps,cur_heatmap), dim=c(11,11,iter+1)) 
	#image(d)
	iter <- iter+1
}
recom_array <- array(c(recom_array,allheatmaps),dim=c(11,11,9,r_iter+1))

}
 
axis_names <- c("2500","5000","10000","20000","40000","80000","160000","320000","640000","1280000","2560000")
dimnames(allheatmaps) <- list(axis_names[11:1], axis_names)

#par(mfrow=c(8,5))
#par(mar=c(5,5,5,5))
#for(i in seq(2,9)){
#	image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
#	axis( 2, at=seq(0,1,length.out=nrow( allheatmaps[,,i] ) ), labels= rownames( allheatmaps[,,i] ), las= 2 )
#	axis( 1, at=seq(0,1,length.out=ncol( allheatmaps[,,i] ) ), labels= colnames( allheatmaps[,,i] ), las= 2 )
#}

#par(mar=c(0,0,0,0))
#par(mfrow=c(8,7))
#for(i in seq(2,9)){
#	for(j in seq(2,(r_iter+1))){
#		#image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
#		#image(allheatmaps[,,i,j],col=pal(50), xaxt="n", yaxt="n")
#		image(recom_array[,,i,j],col=pal(50), xaxt="n", yaxt="n",zlim=c(min(recom_array,na.rm=TRUE),max(recom_array,na.rm=TRUE)))
#		box()
#	}
#}

#dev.print(pdf,paste('colsandrows_all', sep=""))


######colors scaled per row
par(mar=c(0,0,0,0))
par(mfrow=c(8,7))
for(i in seq(2,9)){
	for(j in seq(2,(r_iter+1))){
		#image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
		#image(allheatmaps[,,i,j],col=pal(50), xaxt="n", yaxt="n")
		image(recom_array[,,i,j],col=pal(50), xaxt="n", yaxt="n",zlim=c(min(recom_array,na.rm=TRUE),max(recom_array[,,i,],na.rm=TRUE)))
		box()
	}
}

#dev.print(pdf,paste('colsandrows_all_rf_scaled_per_row.pdf', sep=""))



######colors scale to max possible rf dist
#ntips = num_species(3) * samples per species(8 to 1)
#maxRF dist = 2* (ntips(including outgroup) - 3)

par(mar=c(0,0,0,0))
par(mfrow=c(8,7))
for(i in seq(2,9)){
	numTips <- (3*(10-i))
	maxRFdist = 2* (numTips-3)
	print(maxRFdist)
	for(j in seq(2,(r_iter+1))){
		#image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
		#image(allheatmaps[,,i,j],col=pal(50), xaxt="n", yaxt="n")
		print(max(recom_array[,,i,],na.rm=TRUE))
		image(recom_array[,,i,j],col=pal(50), xaxt="n", yaxt="n",zlim=c(min(recom_array,na.rm=TRUE),maxRFdist))
		box()
	}
}

par(mar=c(0,0,0,0))
par(mfrow=c(8,1))
for(i in seq(2,9)){
	numTips <- (3*(10-i))+1
	maxRFdist = 2* (numTips-3)
	print(maxRFdist)
	#for(j in seq(2,(r_iter+1))){
		#image(t(apply(allheatmaps[,,i],2,rev)),xaxt="n", yaxt="n")
		#image(allheatmaps[,,i,j],col=pal(50), xaxt="n", yaxt="n")
		#print(max(recom_array[,,i,],na.rm=TRUE))
		image(recom_array[,,i,3],col=pal(50), xaxt="n", yaxt="n",zlim=c(min(recom_array,na.rm=TRUE),maxRFdist))
		box()
	#}
}
#dev.print(pdf,paste('col_100_resampled_reps_r0.5.pdf', sep=""))


#dev.print(pdf,paste('colsandrows_all_max_poss_RF.pdf', sep=""))

#dev.print(pdf,paste(rrate,'col.pdf', sep=""))

#sqrt_recom_array <- sqrt(recom_array)

par(mar=c(8,8,4,8))
par(cex.lab=1.5)
image(recom_array[,,3,3],axes=F,col='transparent')
#axis(2,at=seq(0,1,l=ncol(recom_array[,,2,5])),labels=axisnames, cex.axis=1.5)
#axis(1,at=seq(0,1,l=ncol(recom_array[,,2,5])),labels=axisnames, cex.axis=1.5)
image.plot(recom_array[,,3,3],col=pal(50),add=T, axis.args=list(cex.axis=1.5), zlim=c(min(recom_array,na.rm=TRUE),max(recom_array,na.rm=TRUE)),legend.args=list(text='Robinson-Foulds Distance', side=4, font=2, line=2.5, cex=1.3))
box()
title(xlab="Divergence Time T1 (thousands of generations)")
title(ylab="Divergence Time T1 (thousands of generations)")

dev.print(pdf,"RF4.pdf")

