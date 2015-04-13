onefact_prep_run_connected <- function(samp.data, name.delim, num.prep.batches, num.columns, exp.name, sname.id, stype.id, sbio.id, sfact1.id, samptype.rand, omics.nm){
	## Function for randomizing when rand.type == "prep and run orders" and checkbox for keeping prep.batches in run.orders is TRUE
	# samp.data is the data entered by the user
	# name.delim is the delimiter that should be used to split the names
	# num.prep.batches is the number of prep batches needed
	# num.columns is the number of columns on the machine
	# exp.name is a character vector of the experiment name
	# sname.id is the column of the data where sample names are stored
	# stype.id is the part of the string that indicates sample type
	# sbio.id is the part of the string that indicates biological rep
	# sfact1.id is the part of the string that indicates factor 1
	# samptype.rand is the type of sample to be considered
	# omics.nm is a T/F variable that indicates whether OMICS should be added to the front of new sample names

#### Collecting Initial Information ####	
	# split the sample names by the delimiter specified #
	samp.splits = strsplit(as.character(samp.data[,sname.id]), name.delim)

	# pull sample type information #
	stype = as.character(unlist(lapply(samp.splits, function(x) x[stype.id])))

	# pull bio_rep information #
	sbio = as.character(unlist(lapply(samp.splits, function(x) x[sbio.id])))

	# pull factor 1 information #
	sfact1 = as.character(unlist(lapply(samp.splits, function(x) x[sfact1.id])))

	# make a data.frame with all of this relevant info #
	temp = data.frame(samp.data, Samp_Type = stype, Bio_Rep = sbio, Factor1 = as.character(sfact1))

	# pull only rows with the correct sample type to be randomized #
	data4prep = temp[which(temp$Samp_Type==samptype.rand),]
	
#### Prep Order Randomization ####

	## intial sampling info ##
	# calculate the number of whole samples per factor and batch #
	num.per.batch.whole = table(data4prep$Factor1) %/% num.prep.batches

	# calculate the number of leftover samples per factor #
	num.per.batch.part = table(data4prep$Factor1) %% num.prep.batches
	
	
	# some temporary placeholders #
	samps.allocated = list()
	samps.unallocated = list()
	prep.batch.id = list()
	
	# the number of factor1 levels #
	num.fact1 = length(names(num.per.batch.whole)) 

	# assign samples for the number of samples that can be equally distributed for each factor level #
	for(i in 1:num.fact1){
		temp.f1 = data4prep[which(data4prep$Factor1 == names(num.per.batch.whole)[i]),]
		
		# if both partial and whole samples are greater than zero #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] > 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = temp.f1[-samp.id,]
		}
		# if there are no partial samples #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] == 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = NULL
		}
		# if there are no whole samples #
		if(num.per.batch.whole[i] == 0 & num.per.batch.part[i] > 0){
			prep.batch.id[[i]] = NULL
			samps.allocated[[i]] = NULL
			samps.unallocated[[i]] = temp.f1
		}
	}

	# deal with unallocated samples #
	# if both partial and whole samples are greater than zero #
	if(sum(num.per.batch.part > 0) > 0 & sum(num.per.batch.whole > 0) > 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}
	# if only whole batches #
	if(sum(num.per.batch.whole > 0) > 0 & sum(num.per.batch.part > 0) == 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = NULL
	}
	# if only partial batches #
	if(sum(num.per.batch.whole > 0) == 0 & sum(num.per.batch.part > 0) > 0){
		samp.all = NULL
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}	


	# combine allocated and formerly unallocated #
	samp.data.wbatch = rbind(samp.all, samp.unall)
	
	# add a prep order by batch #
	temp.prep.order = list()
	for(j in 1:num.prep.batches){
		temp.f1.prep = samp.data.wbatch[which(samp.data.wbatch$Prep_Batch == j),]
		temp.prep.order[[j]] = data.frame(temp.f1.prep, Prep_Order = sample(1:nrow(temp.f1.prep), nrow(temp.f1.prep)))
	}
	
	# recombine the data #
	samp.data.wbord = do.call(rbind, temp.prep.order)

	# order the data by batch and prep order #
	samp.data.wbord = samp.data.wbord[order(samp.data.wbord$Prep_Batch, samp.data.wbord$Prep_Order), ]

	# add overall prep order #
	samp.data.wbord = data.frame(samp.data.wbord, Overall_Prep_Order = 1:nrow(samp.data.wbord))
	
	# make a summary of treatment by batch #
	prep.tab = table(samp.data.wbord$Factor1, samp.data.wbord$Prep_Batch)
	
#### Run Order Randomization ####

temp.samp.all = list()
temp.samp.una = list()
for(i in 1:num.prep.batches){
	# pull samples that are in current prep batch
	temp.batch = samp.data.wbord[which(samp.data.wbord$Prep_Batch==i),]
	
	# calculate number of whole samples across columns #
	num.whole.cols = table(temp.batch$Factor1)%/% num.columns
	
	# calculate remaining number of samples #
	num.part.cols = table(temp.batch$Factor1)%% num.columns
	
	run.col.ind = list()
	samp.all = list()
	samp.una = list()
	for(j in 1:length(num.whole.cols)){
		# pull current factor 1 level #
		temp.fact1 = temp.batch[which(temp.batch$Factor1 == names(num.whole.cols)[j]),]
		
		# if there are partial and whole samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] > 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = temp.fact1[-sam.id,]
		}
		# if there are no partial samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] == 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = NULL
		}
		# if there are no whole samples #
		if(num.whole.cols[j] == 0 & num.part.cols[j] > 0){
			run.col.ind[[j]] = NULL
			samp.all[[j]] = NULL
			samp.una[[j]] = temp.fact1
		}
	}

	# deal with unallocated samples #
	# if there are some partial and some whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Run_Batch = i, Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Run_Batch = i, Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	#if there are only whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) == 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Run_Batch = i, Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = NULL
	}
	# if there are only partial samples #
	if(sum(num.whole.cols > 0) == 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = NULL
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Run_Batch = i, Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	
}	

	# combine allocated and formerly unallocated samples #
	temp.col.res = rbind(do.call(rbind, temp.samp.all), do.call(rbind, temp.samp.una))
	
	# assign run order within batch and column #
	dat.runord = list()
	for(i in 1:num.prep.batches){
		dat.temp.ord = list()
		for(j in 1:num.columns){
			t.dat = temp.col.res[which(temp.col.res$Run_Batch ==i & temp.col.res$Block == j),]
			runord = sample(1:nrow(t.dat), nrow(t.dat))
			dat.temp.ord[[j]] = data.frame(t.dat, Run_Order = runord)
		}
		dat.runord[[i]] = do.call(rbind, dat.temp.ord)
	}

	# combine across batches #
	near.final.dat = do.call(rbind, dat.runord)
	# order by batch and column #
	near.final.dat = near.final.dat[order(near.final.dat$Run_Batch, near.final.dat$Block),]
	# assign overall run order #
	near.final.dat = data.frame(near.final.dat, Overall_Run_Order = 1:nrow(near.final.dat))
	
	# create temporary sample names vector #
	# get rid of any extraneous symbols #
	near.final.dat$Factor1 = gsub("[^[:alnum:]]","", near.final.dat$Factor1)
	near.final.dat$Samp_Type = gsub("[^[:alnum:]]","", near.final.dat$Samp_Type)
	near.final.dat$Bio_Rep = gsub("[^[:alnum:]]","", near.final.dat$Bio_Rep)

	temp.names = paste(exp.name, near.final.dat$Factor1, near.final.dat$Bio_Rep, near.final.dat$Samp_Type, sprintf("%03d", near.final.dat$Overall_Prep_Order), sep = "_")
	
	if(omics.nm == TRUE){
		temp.names = paste("OMICS",temp.names, sep = "_")
	}

#### Assemble Results ####
rm.cols = c(grep("Samp_Type",names(near.final.dat)), grep("Bio_Rep", names(near.final.dat)), grep("Factor1", names(near.final.dat)))

final.dat = near.final.dat[,-rm.cols]
final.dat[,sname.id] = temp.names
final.dat = final.dat[order(final.dat$Overall_Prep_Order),]

run.tab = data.frame(table(near.final.dat$Run_Batch, near.final.dat$Block, near.final.dat$Factor1))
names(run.tab) = c("Run_Batch","Block","Factor1","Freq")

return(list(rand.data = final.dat, prep.tab = prep.tab, run.tab = run.tab))
}

######################### NOT CONNECTED ##############################

onefact_prep_run_notconnected <- function(samp.data, name.delim, num.prep.batches, num.run.batches, num.columns, exp.name, sname.id, stype.id, sbio.id, sfact1.id, samptype.rand, omics.nm){
	## Function for randomizing when rand.type == "prep and run orders" and checkbox for keeping prep.batches in run.orders is TRUE
	# samp.data is the data entered by the user
	# name.delim is the delimiter that should be used to split the names
	# num.prep.batches is the number of prep batches needed
	# num.run.batches is the number of run batches needed
	# num.columns is the number of columns on the machine
	# exp.name is a character vector of the experiment name
	# sname.id is the column of the data where sample names are stored
	# stype.id is the part of the string that indicates sample type
	# sbio.id is the part of the string that indicates biological rep
	# sfact1.id is the part of the string that indicates factor 1
	# samptype.rand is the type of sample to be considered
	# omics.nm is a T/F variable that indicates whether OMICS should be added to the front of new sample names

#### Collecting Initial Information ####	
	# split the sample names by the delimiter specified #
	samp.splits = strsplit(as.character(samp.data[,sname.id]), name.delim)

	# pull sample type information #
	stype = as.character(unlist(lapply(samp.splits, function(x) x[stype.id])))

	# pull bio_rep information #
	sbio = as.character(unlist(lapply(samp.splits, function(x) x[sbio.id])))

	# pull factor 1 information #
	sfact1 = as.character(unlist(lapply(samp.splits, function(x) x[sfact1.id])))

	# make a data.frame with all of this relevant info #
	temp = data.frame(samp.data, Samp_Type = stype, Bio_Rep = sbio, Factor1 = as.character(sfact1))

	# pull only rows with the correct sample type to be randomized #
	data4prep = temp[which(temp$Samp_Type==samptype.rand),]
	
#### Prep Order Randomization ####

	## intial sampling info ##
	# calculate the number of whole samples per factor and batch #
	num.per.batch.whole = table(data4prep$Factor1) %/% num.prep.batches

	# calculate the number of leftover samples per factor #
	num.per.batch.part = table(data4prep$Factor1) %% num.prep.batches
	
	
	# some temporary placeholders #
	samps.allocated = list()
	samps.unallocated = list()
	prep.batch.id = list()
	
	# the number of factor1 levels #
	num.fact1 = length(names(num.per.batch.whole)) 

	# assign samples for the number of samples that can be equally distributed for each factor level #
	for(i in 1:num.fact1){
		temp.f1 = data4prep[which(data4prep$Factor1 == names(num.per.batch.whole)[i]),]
		
		# if both partial and whole samples are greater than zero #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] > 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = temp.f1[-samp.id,]
		}
		# if there are no partial samples #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] == 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = NULL
		}
		# if there are no whole samples #
		if(num.per.batch.whole[i] == 0 & num.per.batch.part[i] > 0){
			prep.batch.id[[i]] = NULL
			samps.allocated[[i]] = NULL
			samps.unallocated[[i]] = temp.f1
		}
	}

	# deal with unallocated samples #
	# if both partial and whole samples are greater than zero #
	if(sum(num.per.batch.part > 0) > 0 & sum(num.per.batch.whole > 0) > 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}
	# if only whole batches #
	if(sum(num.per.batch.whole > 0) > 0 & sum(num.per.batch.part > 0) == 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = NULL
	}
	# if only partial batches #
	if(sum(num.per.batch.whole > 0) == 0 & sum(num.per.batch.part > 0) > 0){
		samp.all = NULL
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}	


	# combine allocated and formerly unallocated #
	samp.data.wbatch = rbind(samp.all, samp.unall)
	
	# add a prep order by batch #
	temp.prep.order = list()
	for(j in 1:num.prep.batches){
		temp.f1.prep = samp.data.wbatch[which(samp.data.wbatch$Prep_Batch == j),]
		temp.prep.order[[j]] = data.frame(temp.f1.prep, Prep_Order = sample(1:nrow(temp.f1.prep), nrow(temp.f1.prep)))
	}
	
	# recombine the data #
	samp.data.wbord = do.call(rbind, temp.prep.order)

	# order the data by batch and prep order #
	samp.data.wbord = samp.data.wbord[order(samp.data.wbord$Prep_Batch, samp.data.wbord$Prep_Order), ]

	# add overall prep order #
	samp.data.wbord = data.frame(samp.data.wbord, Overall_Prep_Order = 1:nrow(samp.data.wbord))
	
	# make a summary of treatment by batch #
	prep.tab = table(samp.data.wbord$Factor1, samp.data.wbord$Prep_Batch)
	
#### Run Order Randomization ####

num.whole.bats = table(samp.data.wbord$Factor1)%/% num.run.batches
num.part.bats = table(samp.data.wbord$Factor1)%% num.run.batches

temp.batch.all = list()
temp.batch.una = list()
temp.id = list()
for(i in 1:length(num.whole.bats)){
	# pull samples that are in current Factor1 group #
	temp.fact1 = samp.data.wbord[which(samp.data.wbord$Factor1==names(num.whole.bats)[i]),]
	
	if(num.whole.bats[i] > 0 & num.part.bats[i] > 0){
		samps = sample(1:nrow(temp.fact1), num.run.batches*num.whole.bats[i])
		temp.id[[i]] = rep(1:num.run.batches, each = num.whole.bats[i])
		temp.batch.all[[i]] = temp.fact1[samps,]
		temp.batch.una[[i]] = temp.fact1[-samps,]
	}
	if(num.whole.bats[i] > 0 & num.part.bats[i] == 0){
		samps = sample(1:nrow(temp.fact1), num.run.batches*num.whole.bats[i])
		temp.id[[i]] = rep(1:num.run.batches, each = num.whole.bats[i])
		temp.batch.all[[i]] = temp.fact1[samps,]
		temp.batch.una[[i]] = NULL
	}
	if(num.whole.bats[i] == 0 & num.part.bats[i] > 0){
		temp.id[[i]] = NULL
		temp.batch.all[[i]] = NULL
		temp.batch.una[[i]] = temp.fact1[-samps,]
	}
}

	# deal with unallocated samples #
	# if there are some partial and some whole samples #
	if(sum(num.whole.bats > 0) > 0 & sum(num.part.bats > 0) > 0){
		temp.wrun.bat1 = data.frame(do.call(rbind, temp.batch.all), Run_Batch = unlist(temp.id), row.names = NULL)
		temp.wrun.bat2 = data.frame(do.call(rbind, temp.batch.una), Run_Batch = rep(1:num.run.batches, length = nrow(do.call(rbind, temp.batch.una))), row.names = NULL)
		temp.wrun.bat = rbind(temp.wrun.bat1, temp.wrun.bat2)
	}
	#if there are only whole samples #
	if(sum(num.whole.bats > 0) > 0 & sum(num.part.bats > 0) == 0){
		temp.wrun.bat = data.frame(do.call(rbind, temp.batch.all), Run_Batch = unlist(temp.id), row.names = NULL)
	}
	# if there are only partial samples #
	if(sum(num.whole.bats > 0) == 0 & sum(num.part.bats > 0) > 0){
		temp.wrun.bat = data.frame(do.call(rbind, temp.batch.una), Run_Batch = rep(1:num.run.batches, length = nrow(do.call(rbind, temp.batch.una))), row.names = NULL)
	}
	
temp.samp.all = list()
temp.samp.una = list()
	for(i in 1:num.run.batches){
		# pull samples that are in current prep batch
		temp.batch = temp.wrun.bat[which(temp.wrun.bat$Run_Batch==i),]

		# calculate number of whole samples across columns #
		num.whole.cols = table(temp.batch$Factor1)%/% num.columns

		# calculate remaining number of samples #
		num.part.cols = table(temp.batch$Factor1)%% num.columns	
	
	run.col.ind = list()
	samp.all = list()
	samp.una = list()
	for(j in 1:length(num.whole.cols)){
		# pull current factor 1 level #
		temp.fact1 = temp.batch[which(temp.batch$Factor1 == names(num.whole.cols)[j]),]
		
		# if there are partial and whole samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] > 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = temp.fact1[-sam.id,]
		}
		# if there are no partial samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] == 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = NULL
		}
		# if there are no whole samples #
		if(num.whole.cols[j] == 0 & num.part.cols[j] > 0){
			run.col.ind[[j]] = NULL
			samp.all[[j]] = NULL
			samp.una[[j]] = temp.fact1
		}
	}

	# deal with unallocated samples #
	# if there are some partial and some whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	#if there are only whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) == 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = NULL
	}
	# if there are only partial samples #
	if(sum(num.whole.cols > 0) == 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = NULL
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	
}	

	# combine allocated and formerly unallocated samples #
	temp.col.res = rbind(do.call(rbind, temp.samp.all), do.call(rbind, temp.samp.una))
	
	# assign run order within batch and column #
	dat.runord = list()
	for(i in 1:num.run.batches){
		dat.temp.ord = list()
		for(j in 1:num.columns){
			t.dat = temp.col.res[which(temp.col.res$Run_Batch ==i & temp.col.res$Block == j),]
			runord = sample(1:nrow(t.dat), nrow(t.dat))
			dat.temp.ord[[j]] = data.frame(t.dat, Run_Order = runord)
		}
		dat.runord[[i]] = do.call(rbind, dat.temp.ord)
	}

	# combine across batches #
	near.final.dat = do.call(rbind, dat.runord)
	# order by batch and column #
	near.final.dat = near.final.dat[order(near.final.dat$Run_Batch, near.final.dat$Block),]
	# assign overall run order #
	near.final.dat = data.frame(near.final.dat, Overall_Run_Order = 1:nrow(near.final.dat))
	
	# create temporary sample names vector #
	# get rid of any extraneous symbols #
	near.final.dat$Factor1 = gsub("[^[:alnum:]]","", near.final.dat$Factor1)
	near.final.dat$Samp_Type = gsub("[^[:alnum:]]","", near.final.dat$Samp_Type)
	near.final.dat$Bio_Rep = gsub("[^[:alnum:]]","", near.final.dat$Bio_Rep)

	temp.names = paste(exp.name, near.final.dat$Factor1, near.final.dat$Bio_Rep, near.final.dat$Samp_Type, sprintf("%03d", near.final.dat$Overall_Prep_Order), sep = "_")
	
	if(omics.nm == TRUE){
		temp.names = paste("OMICS",temp.names, sep = "_")
	}

#### Assemble Results ####
rm.cols = c(grep("Samp_Type",names(near.final.dat)), grep("Bio_Rep", names(near.final.dat)), grep("Factor1", names(near.final.dat)))

final.dat = near.final.dat[,-rm.cols]
final.dat[,sname.id] = temp.names
final.dat = final.dat[order(final.dat$Overall_Prep_Order),]

run.tab = data.frame(table(near.final.dat$Run_Batch, near.final.dat$Block, near.final.dat$Factor1))
names(run.tab) = c("Run_Batch","Block","Factor1","Freq")

return(list(rand.data = final.dat, prep.tab = prep.tab, run.tab = run.tab))
}

############################# Prep Order Only #############################
onefact_prep <- function(samp.data, name.delim, num.prep.batches, exp.name, sname.id, stype.id, sbio.id, sfact1.id, samptype.rand, omics.nm){
	## Function for randomizing when rand.type == "prep and run orders" and checkbox for keeping prep.batches in run.orders is TRUE
	# samp.data is the data entered by the user
	# name.delim is the delimiter that should be used to split the names
	# num.prep.batches is the number of prep batches needed
	# num.columns is the number of columns on the machine
	# exp.name is a character vector of the experiment name
	# sname.id is the column of the data where sample names are stored
	# stype.id is the part of the string that indicates sample type
	# sbio.id is the part of the string that indicates biological rep
	# sfact1.id is the part of the string that indicates factor 1
	# samptype.rand is the type of sample to be considered
	# omics.nm is a T/F variable that indicates whether OMICS should be added to the front of new sample names

#### Collecting Initial Information ####	
	# split the sample names by the delimiter specified #
	samp.splits = strsplit(as.character(samp.data[,sname.id]), name.delim)

	# pull sample type information #
	stype = as.character(unlist(lapply(samp.splits, function(x) x[stype.id])))

	# pull bio_rep information #
	sbio = as.character(unlist(lapply(samp.splits, function(x) x[sbio.id])))

	# pull factor 1 information #
	sfact1 = as.character(unlist(lapply(samp.splits, function(x) x[sfact1.id])))

	# make a data.frame with all of this relevant info #
	temp = data.frame(samp.data, Samp_Type = stype, Bio_Rep = sbio, Factor1 = as.character(sfact1))

	# pull only rows with the correct sample type to be randomized #
	data4prep = temp[which(temp$Samp_Type==samptype.rand),]
	
#### Prep Order Randomization ####

	## intial sampling info ##
	# calculate the number of whole samples per factor and batch #
	num.per.batch.whole = table(data4prep$Factor1) %/% num.prep.batches

	# calculate the number of leftover samples per factor #
	num.per.batch.part = table(data4prep$Factor1) %% num.prep.batches
	
	
	# some temporary placeholders #
	samps.allocated = list()
	samps.unallocated = list()
	prep.batch.id = list()
	
	# the number of factor1 levels #
	num.fact1 = length(names(num.per.batch.whole)) 

	# assign samples for the number of samples that can be equally distributed for each factor level #
	for(i in 1:num.fact1){
		temp.f1 = data4prep[which(data4prep$Factor1 == names(num.per.batch.whole)[i]),]
		
		# if both partial and whole samples are greater than zero #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] > 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = temp.f1[-samp.id,]
		}
		# if there are no partial samples #
		if(num.per.batch.whole[i] > 0 & num.per.batch.part[i] == 0){
			samp.id = sample(1:nrow(temp.f1), num.prep.batches*num.per.batch.whole[i])
			prep.batch.id[[i]] = rep(1:num.prep.batches, each = num.per.batch.whole[i])
			samps.allocated[[i]] = temp.f1[samp.id,]
			samps.unallocated[[i]] = NULL
		}
		# if there are no whole samples #
		if(num.per.batch.whole[i] == 0 & num.per.batch.part[i] > 0){
			prep.batch.id[[i]] = NULL
			samps.allocated[[i]] = NULL
			samps.unallocated[[i]] = temp.f1
		}
	}

	# deal with unallocated samples #
	# if both partial and whole samples are greater than zero #
	if(sum(num.per.batch.part > 0) > 0 & sum(num.per.batch.whole > 0) > 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}
	# if only whole batches #
	if(sum(num.per.batch.whole > 0) > 0 & sum(num.per.batch.part > 0) == 0){
		samp.all = data.frame(do.call(rbind, samps.allocated), Prep_Batch = unlist(prep.batch.id))
		samp.unall = NULL
	}
	# if only partial batches #
	if(sum(num.per.batch.whole > 0) == 0 & sum(num.per.batch.part > 0) > 0){
		samp.all = NULL
		samp.unall = data.frame(do.call(rbind, samps.unallocated), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samps.unallocated))))
	}	


	# combine allocated and formerly unallocated #
	samp.data.wbatch = rbind(samp.all, samp.unall)
	
	# add a prep order by batch #
	temp.prep.order = list()
	for(j in 1:num.prep.batches){
		temp.f1.prep = samp.data.wbatch[which(samp.data.wbatch$Prep_Batch == j),]
		temp.prep.order[[j]] = data.frame(temp.f1.prep, Prep_Order = sample(1:nrow(temp.f1.prep), nrow(temp.f1.prep)))
	}
	
	# recombine the data #
	samp.data.wbord = do.call(rbind, temp.prep.order)

	# order the data by batch and prep order #
	samp.data.wbord = samp.data.wbord[order(samp.data.wbord$Prep_Batch, samp.data.wbord$Prep_Order), ]

	# add overall prep order #
	samp.data.wbord = data.frame(samp.data.wbord, Overall_Prep_Order = 1:nrow(samp.data.wbord))
	
	# make a summary of treatment by batch #
	prep.tab = table(samp.data.wbord$Factor1, samp.data.wbord$Prep_Batch)
	
	near.final.dat = samp.data.wbord

		# create temporary sample names vector #
		# get rid of any extraneous symbols #
		near.final.dat$Factor1 = gsub("[^[:alnum:]]","", near.final.dat$Factor1)
		near.final.dat$Samp_Type = gsub("[^[:alnum:]]","", near.final.dat$Samp_Type)
		near.final.dat$Bio_Rep = gsub("[^[:alnum:]]","", near.final.dat$Bio_Rep)

		temp.names = paste(exp.name, near.final.dat$Factor1, near.final.dat$Bio_Rep, near.final.dat$Samp_Type, sprintf("%03d", near.final.dat$Overall_Prep_Order), sep = "_")

		if(omics.nm == TRUE){
			temp.names = paste("OMICS",temp.names, sep = "_")
		}

	#### Assemble Results ####
	rm.cols = c(grep("Samp_Type",names(near.final.dat)), grep("Bio_Rep", names(near.final.dat)), grep("Factor1", names(near.final.dat)))

	final.dat = near.final.dat[,-rm.cols]
	final.dat[,sname.id] = temp.names
	final.dat = final.dat[order(final.dat$Overall_Prep_Order),]

	return(list(rand.data = final.dat, prep.tab = prep.tab))
}

########################## Run Order Only ########################
onefact_run <- function(samp.data, name.delim, num.run.batches, num.columns, exp.name, sname.id, stype.id, sbio.id, sfact1.id, samptype.rand, omics.nm){
	## Function for randomizing when rand.type == "prep and run orders" and checkbox for keeping prep.batches in run.orders is TRUE
	# samp.data is the data entered by the user
	# name.delim is the delimiter that should be used to split the names
	# num.prep.batches is the number of prep batches needed
	# num.run.batches is the number of run batches needed
	# num.columns is the number of columns on the machine
	# exp.name is a character vector of the experiment name
	# sname.id is the column of the data where sample names are stored
	# stype.id is the part of the string that indicates sample type
	# sbio.id is the part of the string that indicates biological rep
	# sfact1.id is the part of the string that indicates factor 1
	# samptype.rand is the type of sample to be considered
	# omics.nm is a T/F variable that indicates whether OMICS should be added to the front of new sample names

#### Collecting Initial Information ####	
	# split the sample names by the delimiter specified #
	samp.splits = strsplit(as.character(samp.data[,sname.id]), name.delim)

	# pull sample type information #
	stype = as.character(unlist(lapply(samp.splits, function(x) x[stype.id])))

	# pull bio_rep information #
	sbio = as.character(unlist(lapply(samp.splits, function(x) x[sbio.id])))

	# pull factor 1 information #
	sfact1 = as.character(unlist(lapply(samp.splits, function(x) x[sfact1.id])))

	# make a data.frame with all of this relevant info #
	temp = data.frame(samp.data, Samp_Type = stype, Bio_Rep = sbio, Factor1 = as.character(sfact1))

	# pull only rows with the correct sample type to be randomized #
	data4prep = temp[which(temp$Samp_Type==samptype.rand),]
	
	
#### Run Order Randomization ####

num.whole.bats = table(data4prep$Factor1)%/% num.run.batches
num.part.bats = table(data4prep$Factor1)%% num.run.batches

temp.batch.all = list()
temp.batch.una = list()
temp.id = list()
for(i in 1:length(num.whole.bats)){
	# pull samples that are in current Factor1 group #
	temp.fact1 = data4prep[which(data4prep$Factor1==names(num.whole.bats)[i]),]
	
	if(num.whole.bats[i] > 0 & num.part.bats[i] > 0){
		samps = sample(1:nrow(temp.fact1), num.run.batches*num.whole.bats[i])
		temp.id[[i]] = rep(1:num.run.batches, each = num.whole.bats[i])
		temp.batch.all[[i]] = temp.fact1[samps,]
		temp.batch.una[[i]] = temp.fact1[-samps,]
	}
	if(num.whole.bats[i] > 0 & num.part.bats[i] == 0){
		samps = sample(1:nrow(temp.fact1), num.run.batches*num.whole.bats[i])
		temp.id[[i]] = rep(1:num.run.batches, each = num.whole.bats[i])
		temp.batch.all[[i]] = temp.fact1[samps,]
		temp.batch.una[[i]] = NULL
	}
	if(num.whole.bats[i] == 0 & num.part.bats[i] > 0){
		temp.id[[i]] = NULL
		temp.batch.all[[i]] = NULL
		temp.batch.una[[i]] = temp.fact1[-samps,]
	}
}

	# deal with unallocated samples #
	# if there are some partial and some whole samples #
	if(sum(num.whole.bats > 0) > 0 & sum(num.part.bats > 0) > 0){
		temp.wrun.bat1 = data.frame(do.call(rbind, temp.batch.all), Run_Batch = unlist(temp.id), row.names = NULL)
		temp.wrun.bat2 = data.frame(do.call(rbind, temp.batch.una), Run_Batch = rep(1:num.run.batches, length = nrow(do.call(rbind, temp.batch.una))), row.names = NULL)
		temp.wrun.bat = rbind(temp.wrun.bat1, temp.wrun.bat2)
	}
	#if there are only whole samples #
	if(sum(num.whole.bats > 0) > 0 & sum(num.part.bats > 0) == 0){
		temp.wrun.bat = data.frame(do.call(rbind, temp.batch.all), Run_Batch = unlist(temp.id), row.names = NULL)
	}
	# if there are only partial samples #
	if(sum(num.whole.bats > 0) == 0 & sum(num.part.bats > 0) > 0){
		temp.wrun.bat = data.frame(do.call(rbind, temp.batch.una), Run_Batch = rep(1:num.run.batches, length = nrow(do.call(rbind, temp.batch.una))), row.names = NULL)
	}
	
temp.samp.all = list()
temp.samp.una = list()
	for(i in 1:num.run.batches){
		# pull samples that are in current prep batch
		temp.batch = temp.wrun.bat[which(temp.wrun.bat$Run_Batch==i),]

		# calculate number of whole samples across columns #
		num.whole.cols = table(temp.batch$Factor1)%/% num.columns

		# calculate remaining number of samples #
		num.part.cols = table(temp.batch$Factor1)%% num.columns	
	
	run.col.ind = list()
	samp.all = list()
	samp.una = list()
	for(j in 1:length(num.whole.cols)){
		# pull current factor 1 level #
		temp.fact1 = temp.batch[which(temp.batch$Factor1 == names(num.whole.cols)[j]),]
		
		# if there are partial and whole samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] > 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = temp.fact1[-sam.id,]
		}
		# if there are no partial samples #
		if(num.whole.cols[j] > 0 & num.part.cols[j] == 0){
			sam.id = sample(1:nrow(temp.fact1), num.columns*num.whole.cols[j])
			run.col.ind[[j]] = rep(1:num.columns, each = num.whole.cols[j])
			samp.all[[j]] = temp.fact1[sam.id,]
			samp.una[[j]] = NULL
		}
		# if there are no whole samples #
		if(num.whole.cols[j] == 0 & num.part.cols[j] > 0){
			run.col.ind[[j]] = NULL
			samp.all[[j]] = NULL
			samp.una[[j]] = temp.fact1
		}
	}

	# deal with unallocated samples #
	# if there are some partial and some whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	#if there are only whole samples #
	if(sum(num.whole.cols > 0) > 0 & sum(num.part.cols > 0) == 0){
		temp.samp.all[[i]] = data.frame(do.call(rbind, samp.all), Block = unlist(run.col.ind), row.names = NULL)
		temp.samp.una[[i]] = NULL
	}
	# if there are only partial samples #
	if(sum(num.whole.cols > 0) == 0 & sum(num.part.cols > 0) > 0){
		temp.samp.all[[i]] = NULL
		temp.samp.una[[i]] = data.frame(do.call(rbind, samp.una), Block = rep(1:num.columns, length = nrow(do.call(rbind, samp.una))), row.names = NULL)
	}
	
}	

	# combine allocated and formerly unallocated samples #
	temp.col.res = rbind(do.call(rbind, temp.samp.all), do.call(rbind, temp.samp.una))
	
	# assign run order within batch and column #
	dat.runord = list()
	for(i in 1:num.run.batches){
		dat.temp.ord = list()
		for(j in 1:num.columns){
			t.dat = temp.col.res[which(temp.col.res$Run_Batch ==i & temp.col.res$Block == j),]
			runord = sample(1:nrow(t.dat), nrow(t.dat))
			dat.temp.ord[[j]] = data.frame(t.dat, Run_Order = runord)
		}
		dat.runord[[i]] = do.call(rbind, dat.temp.ord)
	}

	# combine across batches #
	near.final.dat = do.call(rbind, dat.runord)
	# order by batch and column #
	near.final.dat = near.final.dat[order(near.final.dat$Run_Batch, near.final.dat$Block, near.final.dat$Run_Order),]
	# assign overall run order #
	near.final.dat = data.frame(near.final.dat, Overall_Run_Order = 1:nrow(near.final.dat))
	
	# create temporary sample names vector #
	# get rid of any extraneous symbols #
	near.final.dat$Factor1 = gsub("[^[:alnum:]]","", near.final.dat$Factor1)
	near.final.dat$Samp_Type = gsub("[^[:alnum:]]","", near.final.dat$Samp_Type)
	near.final.dat$Bio_Rep = gsub("[^[:alnum:]]","", near.final.dat$Bio_Rep)

	temp.names = paste(exp.name, near.final.dat$Factor1, near.final.dat$Bio_Rep, near.final.dat$Samp_Type, sprintf("%03d", near.final.dat$Overall_Run_Order), sep = "_")
	
	if(omics.nm == TRUE){
		temp.names = paste("OMICS",temp.names, sep = "_")
	}

#### Assemble Results ####
rm.cols = c(grep("Samp_Type",names(near.final.dat)), grep("Bio_Rep", names(near.final.dat)), grep("Factor1", names(near.final.dat)))

final.dat = near.final.dat[,-rm.cols]
final.dat[,sname.id] = temp.names
final.dat = final.dat[order(final.dat$Overall_Run_Order),]

run.tab = data.frame(table(near.final.dat$Run_Batch, near.final.dat$Block, near.final.dat$Factor1))
names(run.tab) = c("Run_Batch","Block","Factor1","Freq")

return(list(rand.data = final.dat, run.tab = run.tab))
}