#Author Christopher Rodriguez

file_extension = "fas"	#File ending (the "fas" in e1_81.fas)
tolerance_list = c("0.01")	#The difference tolerances you want to run the hivnet_cluster at
do_mafft_align = FALSE	#Do you want to run an alignment of the sequences? 
should_remove_nonclustered_seqs = FALSE	#Should we remove the sequences that don't cluster from the fasta file? 
should_make_stats = TRUE #ld we produce the stat csv file that describes our clustering?
should_make_roc_curves = TRUE	#Should we make roc curves of each fasta?
if(should_make_roc_curves)
	library("pROC")
	
patientIDLoc = 1	#Where the patientID is located: 
					#e.g. in 10014_8_1a_20031029, 10014 is
					#the patientID, so patientIDLoc is set to 1.
					#Used for roc curves (inter v intra patients)



#Writes a statistical summary of the data obtained from the cluster run
#Function contains code written by James Jarad Dollar
write_cluster_stats <- function(f_name)
{
	cluster_data = read.csv(file = paste(f_name,"_cluster.csv",sep=""))
	stat_data = data.frame(row.names = c("Total Number of Sequences","Number in Clusters", "Number of Clusters", "Average Number in Clusters", "Maximum Number in Clusters", "Minimum Number in Clusters", "1st Quartile", "3rd Quartile", "Number of Monophyletic Clusters", "Percent Monophyletic Clusters",  "Number of Non-Monophyletic Clusters", "Percent Non-Monophyletic Clusters"))
	for(i in 2:ncol(cluster_data))
	{
		clusters <- table(cluster_data[which(cluster_data[,i]!=0),i])
		nclust <- c(length(clusters))
		aveclust <- c(mean(clusters))
		maxclust <- c(max(clusters)) 
		minclust <- c(min(clusters)) 
		quart1 <- as.numeric(quantile(clusters, .25))
		quart3 <- as.numeric(quantile(clusters, .75))
		twoclust <- c(length(which(clusters == 2)))
		gtwoclust <- c(length(which(clusters > 2)))
		tperclust <- c((length(which(clusters == 2))/length(clusters))*100)
		gtperclust <- c((length(which(clusters > 2))/length(clusters))*100)
		stat_data = cbind(stat_data,c(sum(clusters)+length(which(cluster_data[,i]==0)),sum(clusters), nclust, aveclust, maxclust, minclust, quart1, quart3, twoclust, tperclust, gtwoclust, gtperclust))
		colnames(stat_data)[i-1] = colnames(cluster_data)[i]
	}
	write.csv(stat_data,paste(f_name,"_stats.csv",sep=""),quote=FALSE)
}


#Removes a non-clustered sequence from the plotn_fname.fsa file
#Note that this is not the "original" fasta file, but rather the alignment/copy 
remove_sequence_from_original <- function(seq_name,f_name_stub)
{
	#The line below returns the name of the original fasta sequence instead of the plotn_fname.fsa file
	#f_name_original = regmatches(f_name_stub, regexpr("_", f_name_stub), invert = TRUE)[[1]][2]

	f_name = paste(f_name_stub,".",file_extension,sep="")
	print(seq_name)
	print(f_name)
	conn <- file(f_name,open="r")
	linn <-readLines(conn)
	start_seq = -1
	for (i in 1:length(linn))
	{
		if(start_seq!=-1)
		{
			if(grepl(">",linn[i]))
			{
				new_f = linn[-(start_seq:(i-1))]
				close(conn)
				writeLines(new_f,f_name)
				print(i)
				break
			}
		}
		else if(grepl(seq_name,linn[i]))
		{
			start_seq = i
		}
		
	}
}

#Main function for running through files
individual <- function()
{
	files = list.files(pattern = paste("\\.",file_extension,"$",sep=""))
	all_vals = data.frame()
	#remove_fsa_indexes(files)
	write("",file="noncluster_list.dat")
	files = files[which(!grepl("plot",regmatches(files, regexpr("_", files), invert = TRUE)))]
	count = 1
	for(f in files)
	{
		f_name = strsplit(f, paste("\\.",file_extension,sep=""))[[1]][1]
		print(f_name)
		new_name = paste("plot",count,"_",f_name,sep="")
		if(do_mafft_align)
		{
			align_fasta(f_name,count)
		}
		else
		{
			
			file.copy(paste(f_name,".",file_extension,sep=""),paste(new_name,".",file_extension,sep=""),overwrite=TRUE)
		}
		tn93_dist(paste(new_name,".",file_extension,sep=""),new_name)
		if(should_make_roc_curves)
		{
			interpats = roc_plot(new_name)
		}
		for(i in 1:length(tolerance_list))
		{
			first_run = (i==1)
			hivnet_cluster(new_name,tolerance_list[i],first_run)
		}
		if(should_make_stats)
		{
			write_cluster_stats(new_name)
		}
		#makeplot(f_name)
		csv_vals = read.csv(file=paste(new_name,".csv",sep=""))
		if(should_make_roc_curves)
		{
			csv_vals = csv_vals[interpats,]
		}
		csv_dists = csv_vals$Distance
		this_vals = cbind(Fasta = rep(new_name,times = length(csv_dists)),Distance = csv_dists,Index = rep(count,times = length(csv_dists)))
		all_vals = rbind(all_vals,this_vals)
		count = count + 1
	}
	
	write.csv(all_vals,"super_plot.csv",row.names = FALSE)
	make_super_plot()

}


#Plots the ROC curve for intra/inter patient data
#Returns the indexes of inter-patient sequences (not intra-patient)
#1 is inter-patient, 0 is intra-patient
roc_plot<-function(f_name)
{

	csv_dat = read.csv(paste(f_name,".csv",sep=""),stringsAsFactors = FALSE)
	
	responses =c(rep(0,nrow(csv_dat)))
	pat_names = character(nrow(csv_dat))
	pat_dist = c(rep(-1,nrow(csv_dat)))
	count = 1
	for(i in 1:nrow(csv_dat)){
		if(unlist(strsplit(csv_dat$ID1[i],"_"))[patientIDLoc] != unlist(strsplit(csv_dat$ID2[i],"_"))[patientIDLoc])
		{	
			responses[i]=1
		}
		else
		{
			pat_names[count] = unlist(strsplit(csv_dat$ID1[i],"_"))[patientIDLoc]
			pat_dist[count] = csv_dat$Distance[i]
			count = count + 1
		}
	}
	pat_names = pat_names[which(pat_names!="")]
	pat_dist = pat_dist[which(pat_dist!=-1)]
	
	inter_intra_dat =c(rep("Intra",nrow(csv_dat)))
	inter_intra_dat[which(responses==1)] = "Inter"
	
	patdat = data.frame(Patient = pat_names,Distance = pat_dist)
	makeplot_patients(paste(f_name,"_pat_boxplots",sep=""),patdat)
	
	inter_intra_df = data.frame(Patient = inter_intra_dat,Distance = csv_dat$Distance)
	makeplot_inter_intra(paste(f_name,"_inter_intra_boxplots",sep=""),inter_intra_df)

	make_histogram_intra_inter(paste(f_name,"_inter_hist",sep=""),csv_dat$Distance[which(responses==1)],"Inter-Patient Sequences")
	make_histogram_intra_inter(paste(f_name,"_intra_hist",sep=""),csv_dat$Distance[which(responses==0)],"Intra-Patient Sequences")
	
	
	roctest = roc(responses,csv_dat$Distance,levels=c(0,1), direction="<")
	pdf(paste(f_name,"_roc.pdf",sep=""))
	plot(roctest)
	text(0.3,0.3,paste("Best Threshold: ",unname(coords(roctest,"best"))[1],sep=""))
	text(0.3,0.25,paste("AUC: ",auc(roctest),sep=""))
	dev.off()
	
	
	return (which(responses==1))
}

#Renames the plotn_fname.fsa file to just a fname.fsa file
#Not in use right now
remove_fsa_indexes<-function(files)
{
	for(f_name in files)
	{
		if(grepl("plot",regmatches(f_name, regexpr("_", f_name), invert = TRUE)[[1]][1]))
		{
			file.rename(f_name,regmatches(f_name, regexpr("_", f_name), invert = TRUE)[[1]][2])
		}
	}
}

#Runs the mafft command to align a sequence
align_fasta <- function(f_name,count)
{
	new_name = paste("plot",count,"_",f_name,".",file_extension,sep="")
	print(paste("mafft ",f_name,'.',file_extension," > ",new_name,sep=""))
	try(system(paste("mafft ",f_name,'.',file_extension," > ",new_name,sep=""), intern = TRUE))
}


#Runs the tn93 command on a fasta file to get pairwise distances
#infile: give it with the .fsa or .fasta extension
#outfile: give it without the .csv extension  
tn93_dist <- function(infile, outfile)
{
	try(system(paste("tn93 -t 1.0 -o ",outfile,".csv"," ",infile,sep=""), intern = TRUE))
}


#Runs the hivnetworkcsv command with the given tolerance to produce sequence clusters
hivnet_cluster <- function(f_name,tolerance,first_run)
{
	if(first_run)
	{
		temp_ext = ""
	}
	else
	{
		temp_ext = "_tmp"
	}
	
	print(paste("hivnetworkcsv -i ",f_name,".csv"," -f plain -c ",f_name,"_cluster",temp_ext,".csv"," -t ",tolerance,sep=""))
	try(system(paste("hivnetworkcsv -i ",f_name,".csv"," -f plain -c ",f_name,"_cluster",temp_ext,".csv"," -t ",tolerance,sep=""), intern = TRUE))
	
	
	cluster = read.csv(file=paste(f_name,"_cluster",temp_ext,".csv",sep=""))
	normal = read.csv(paste(f_name,".csv",sep=""),stringsAsFactors=FALSE)
	all_seqs = unique(c(normal$ID1,normal$ID2))
	cluster_seqs = unique(cluster$SequenceID)
	differences = unlist(setdiff(all_seqs,cluster_seqs))
	diff_combo = cbind(differences,rep("0",length(differences)))
	colnames(diff_combo) = colnames(cluster)
	ncluster = rbind(cluster,diff_combo)
	ncluster$SequenceID = as.character(ncluster$SequenceID)
	ncluster = ncluster[order(ncluster$SequenceID),]
	if(first_run)
	{
		colnames(ncluster)[length(ncluster)] = paste("Cluster_",tolerance,sep="")
		write.csv(ncluster,paste(f_name,"_cluster.csv",sep=""), row.names = FALSE,quote=FALSE)
	}
	else
	{
		ocluster = read.csv(file=paste(f_name,"_cluster.csv",sep=""))
		ocluster = cbind(ocluster,ncluster$ClusterID)
		colnames(ocluster)[length(ocluster)] = paste("Cluster_",tolerance,sep="")
		write.csv(ocluster,paste(f_name,"_cluster.csv",sep=""), row.names = FALSE,quote=FALSE)
		file.remove(paste(f_name,"_cluster",temp_ext,".csv",sep=""))
	}
	for(d in differences)
	{
		if(should_remove_nonclustered_seqs)
		{
			remove_sequence_from_original(d,f_name)
		}
	}
	write(paste(f_name,"_",tolerance,":",paste(differences,collapse=','),sep=""),file="noncluster_list.dat",append=TRUE)
}



#Creates a boxplot of the pairwise distances of sequences in 1 fasta file 
makeplot <- function(csv_file)
{
	read.csv(file=paste(csv_file,".csv",sep="")) -> x
	pdf(file = paste(csv_file,".pdf",sep=""),width=5+length(table(x$ID1)))
	boxplot(Distance ~ ID1, data = x,main = paste("Pairwise distances of",csv_file), xlab = "Sequence",ylab = "Distance",xaxt="n")
	axis(1,at=1:length(table(x$ID1)), labels=1:length(table(x$ID1)))
	dev.off()
}

#Creates a boxplot of all fastas in a directory
#Think of it as a combination of the individual boxplots
#The labels are 1 through n which correspond to the plotn_fname.fsa file
make_super_plot <- function()
{
	read.csv(file="super_plot.csv") -> dat_frame
	print(length(table(dat_frame[,1])))
	pdf(file = "super_plot.pdf",width = 5+length(table(dat_frame[,1])))
	boxplot(Distance ~ Fasta, data = dat_frame,main = "Pairwise distances of all", xlab = "Fasta",ylab = "Distance",xaxt="n")
	axis(1,at=1:length(table(dat_frame[,1])), labels=1:length(table(dat_frame[,1])))
	dev.off()
	
}

#Creates a boxplot with actual sequence names instead of numbers as labels
#Deprecated, but left in for future reference
makeplot.correctLabels <- function(csv_file)
{
	read.csv(file=paste(csv_file,".csv",sep="")) -> x
	pdf(file = paste(csv_file,".pdf",sep=""),width=10+length(table(x$ID1))*4)
	boxplot(Distance ~ ID1, data = x, main = paste("Pairwise distances of",csv_file), xlab = "Sequence",boxwex = 0.4,ylab = "Distance",cex.axis=0.5,yaxt="n")
	axis(2,cex.axis=1)
	dev.off()
}


#Creates a boxplot of the distances for each patient followed by a boxplot for intra-patient distances and a boxplot for inter-patient distances
makeplot_patients <- function(f_name,patdat)
{
	pdf(file = paste(f_name,'.pdf',sep=""),width = 5+length(table(patdat[,1])))
	boxplot(Distance ~ Patient, data = patdat,main = "Distances of Patients", xlab = "Patient",ylab = "Distance")
	dev.off()
}
makeplot_inter_intra <- function(f_name,patdat)
{
	pdf(file = paste(f_name,'.pdf',sep=""),width = 5+length(table(patdat[,1])))
	boxplot(Distance ~ Patient, data = patdat,main = "Distances of Inter-Patient Sequences v. Intra-Patient Sequences", xlab = "Type",ylab = "Distance")
	dev.off()
}

make_histogram_intra_inter <- function(f_name,dat,type)
{
	pdf(file = paste(f_name,'.pdf',sep=""))
	hist(dat,main = paste("Distances of ",type,sep=""), xlab = "Distance",ylab = "Frequency")
	dev.off()
}

individual()