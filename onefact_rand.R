
onefact_runprep_random <- function(samp.data, name.delim, rand.type, pr.link = FALSE, num.prep.batches = 0, num.run.batches=0,num.columns = 0, exp.name, sname.id, stype.id, sbio.id, sfact1.id, samptype.rand,omics.nm){

	samp.splits = strsplit(as.character(samp.data[,sname.id]), name.delim)
	stype = as.character(unlist(lapply(samp.splits, function(x) x[stype.id])))
	sbio = as.character(unlist(lapply(samp.splits, function(x) x[sbio.id])))
	sfact1 = as.character(unlist(lapply(samp.splits, function(x) x[sfact1.id])))
	sfact1 = sfact1[which(stype==samptype.rand)]
	sampdata.wid = data.frame(samp.data[which(stype==samptype.rand),], ID = 1:sum(stype==samptype.rand))

if(rand.type == "prep and run orders"){
		
		###### Prep Order #########
		# calculate batch sizes #
		numpfact.whole = table(sfact1) %/% num.prep.batches
		numpfact.part = table(sfact1) %% num.prep.batches
		
		samp.alloc = list()
		samp.unalloc = list()
		prep.batch.ind = list()
		for(i in 1:length(names(numpfact.whole))){
				temp = sampdata.wid[which(sfact1 == names(numpfact.whole[i])),]
				if(numpfact.whole[i] > 0 & numpfact.part[i] > 0){
				samp.id = sample(1:nrow(temp), num.prep.batches*numpfact.whole[i])
				prep.batch.ind[[i]] = rep(1:num.prep.batches, each = numpfact.whole[i])
				samp.alloc[[i]] = temp[samp.id,]
				samp.unalloc[[i]]= temp[-samp.id,]
				}
				if(numpfact.whole[i] > 0 & numpfact.part[i] == 0){
				samp.id = sample(1:nrow(temp), num.prep.batches*numpfact.whole[i])
				prep.batch.ind[[i]] = rep(1:num.prep.batches, each = numpfact.whole[i])
				samp.alloc[[i]] = temp[samp.id,]
				samp.unalloc[[i]]= NULL
				}
				if(numpfact.whole[i] == 0 & numpfact.part[i] > 0){
				prep.batch.ind[[i]] = NULL
				samp.alloc[[i]] = NULL
				samp.unalloc[[i]] = temp
				}
			}
		
		
		if(sum(numpfact.whole > 0)>0 & sum(numpfact.part > 0)>0){
			samp.all = data.frame(do.call(rbind, samp.alloc), Prep_Batch = unlist(prep.batch.ind))
			samp.unall = cbind(do.call(rbind, samp.unalloc), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samp.unalloc))))
		}
		if(sum(numpfact.whole > 0)>0 & sum(numpfact.part > 0)==0){
			samp.all = data.frame(do.call(rbind, samp.alloc), Prep_Batch = unlist(prep.batch.ind))
			samp.unall = NULL
		}
		if(sum(numpfact.whole > 0)==0 & sum(numpfact.part > 0)>0){
			samp.all = NULL
			samp.unall = cbind(do.call(rbind, samp.unalloc), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samp.unalloc))))
		}

		samp.dat.wbatch = rbind(samp.all,samp.unall)
	
		temp.pord = list()
		for(j in 1:num.prep.batches){
			temp.dat = samp.dat.wbatch[which(samp.dat.wbatch$Prep_Batch==j),]
			temp.pord[[j]] = data.frame(temp.dat, Prep_Order = sample(1:nrow(temp.dat),nrow(temp.dat)))
		}
		
		samp.dat.wbatord = do.call(rbind, temp.pord)
		samp.dat.wbatord = samp.dat.wbatord[order(samp.dat.wbatord$Prep_Batch, samp.dat.wbatord$Prep_Order),]
		samp.dat.wbatord = data.frame(samp.dat.wbatord, Overall_Prep_Order = 1:nrow(samp.dat.wbatord))
		
		junk = data.frame(Factor = sfact1[samp.dat.wbatord$ID], Batch = samp.dat.wbatord$Prep_Batch)
		prep.tab = table(junk$Factor, junk$Batch)
		
		### if prep and run are linked ###
		if(pr.link == TRUE){
			temp.col.res.a = list()
			temp.col.res.u = list()
		for(i in 1:num.prep.batches){
			temp.d = samp.dat.wbatord[which(samp.dat.wbatord$Prep_Batch == i),]
			temp.fact = sfact1[temp.d$ID]
			
			numpf.whole = table(temp.fact) %/% num.columns
			numpf.part = table(temp.fact) %% num.columns

			samp.a = list()
			samp.u = list()
			run.col.ind = list()
			for(j in 1:length(names(numpf.whole))){
					tem = temp.d[which(temp.fact == names(numpf.whole[j])),]
				if(numpf.whole[j] > 0 & numpf.part[j] >0){
					sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
					run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
					samp.a[[j]] = tem[sam.id,]
					samp.u[[j]]= tem[-sam.id,]
				}
				if(numpf.whole[j] == 0 & numpf.part[j] >0){
						samp.a[[j]] = NULL
						samp.u[[j]] = tem
				}
				if(numpf.whole[j] > 0 & numpf.part[j] == 0){
					sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
					run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
					samp.a[[j]] = tem[sam.id,]
					samp.u[[j]]= NULL
				}	
			}
			
			if(sum(numpf.whole > 0) > 0 & sum(numpf.part > 0) > 0){
			temp.col.res.a[[i]] = data.frame(do.call(rbind, samp.a), Run_Batch = i, Column = unlist(run.col.ind))
			temp.col.res.u[[i]] = data.frame(do.call(rbind, samp.u), Run_Batch = i, Column = rep(1:num.columns, length = nrow(do.call(rbind, samp.u))))
			}
			if(sum(numpf.whole > 0) == 0 & sum(numpf.part > 0) > 0){
			temp.col.res.a[[i]] = NULL
			temp.col.res.u[[i]] = data.frame(do.call(rbind, samp.u), Run_Batch = i, Column = rep(1:num.columns, length = nrow(do.call(rbind, samp.u))))
			}
			if(sum(numpf.whole > 0) > 0 & sum(numpf.part > 0) == 0){
			temp.col.res.a[[i]] = data.frame(do.call(rbind, samp.a), Run_Batch = i, Column = unlist(run.col.ind))
			temp.col.res.u[[i]] = NULL
			}

		}
		
		temp.run.res = rbind(do.call(rbind, temp.col.res.a), do.call(rbind, temp.col.res.u))
			
		temp.fina = list()
		for(i in 1:num.prep.batches){
			temp.fin = list()
			for(j in 1:num.columns){
				tt = temp.run.res[which(temp.run.res$Run_Batch == i & temp.run.res$Column==j),]
				runord = sample(1:nrow(tt),nrow(tt))
				temp.fin[[j]] = data.frame(tt, Run_Order = runord)
				}
				temp.fina[[i]] = do.call(rbind, temp.fin)
		}

		temp.final = do.call(rbind, temp.fina)
		temp.final = temp.final[order(temp.final$Run_Batch,temp.final$Column,temp.final$Run_Order),]
		temp.final = data.frame(temp.final, Overall_Run_Order = 1:nrow(temp.final))
		
		junk2 = data.frame(Factor = sfact1[temp.final$ID], Column = temp.final$Column, Batch = temp.final$Run_Batch)
		run.tab = data.frame(table(junk2$Factor, junk2$Column, junk2$Batch))
		
		## assemble new sample names ##
		n.exp = gsub("[^[:alnum:]]","",exp.name)
		n.sfact1 = gsub("[^[:alnum:]]","",sfact1)
		n.stype = gsub("[^[:alnum:]]","",stype)
		n.sbio = gsub("[^[:alnum:]]","",sbio)

		if(omics.nm==TRUE){new.names = paste("OMICS",n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID], sprintf("%03d",temp.final$Overall_Prep_Order), sep="_")}else{
			new.names = paste(n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID],sprintf("%03d",temp.final$Overall_Prep_Order), sep="_")
		}
	
		## assemble results ##
		temp.final[,grep("ID", names(temp.final))] = new.names
		names(temp.final)[grep("ID", names(temp.final))] = "New_Name"
		names(temp.final)[sname.id] = "Original_Name"
			return(list(random = temp.final, tab1 = prep.tab,tab2 = run.tab))
	}else{
				temp.d = samp.dat.wbatord
				temp.fact = sfact1[temp.d$ID]

				numpf.whole = table(temp.fact) %/% num.columns
				numpf.part = table(temp.fact) %% num.columns

				temp.col.res.a = list()
				temp.col.res.u = list()
				run.col.ind = list()
				for(j in 1:length(names(numpf.whole))){
						tem = temp.d[which(temp.fact == names(numpf.whole[j])),]
						if(numpf.whole[j] > 0 & numpf.part[j] > 0){
						sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
						run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
						temp.col.res.a[[j]] = tem[sam.id,]
						temp.col.res.u[[j]]= tem[-sam.id,]
						}
						if(numpf.whole[j] == 0 & numpf.part[j] > 0){
							run.col.ind[[j]] = NULL
							temp.col.res.a[[j]] = NULL
							temp.col.res.u[[j]] = tem
						}
						if(numpf.whole[j] > 0 & numpf.part[j] == 0){
							sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
							run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
							temp.col.res.a[[j]] = tem[sam.id,]
							temp.col.res.u[[j]]= NULL
						}
				}

				if(sum(numpf.whole > 0) > 0 & sum(numpf.part >0) >0){
				temp.col.a = data.frame(do.call(rbind, temp.col.res.a), Column = unlist(run.col.ind))
				temp.col.u = data.frame(do.call(rbind, temp.col.res.u), Column = rep(1:num.columns, length = nrow(do.call(rbind, temp.col.res.u))))
				}
				if(sum(numpf.whole > 0) == 0 & sum(numpf.part >0) >0){
					temp.col.a = NULL
					temp.col.u = data.frame(do.call(rbind, temp.col.res.u), Column = rep(1:num.columns, length = nrow(do.call(rbind, temp.col.res.u))))
				}
				if(sum(numpf.whole > 0) > 0 & sum(numpf.part >0) ==0){
				temp.col.a = data.frame(do.call(rbind, temp.col.res.a), Column = unlist(run.col.ind))
				temp.col.u = NULL
				}
			temp.run.res = rbind(temp.col.a, temp.col.u)

			temp.fina = list()
				for(j in 1:num.columns){
					tt = temp.run.res[which(temp.run.res$Column==j),]
					runord = sample(1:nrow(tt),nrow(tt))
					temp.fina[[j]] = data.frame(tt, Run_Order = runord)
					}

			temp.final = do.call(rbind, temp.fina)
			temp.final = temp.final[order(temp.final$Column,temp.final$Run_Order),]
			temp.final = data.frame(temp.final, Overall_Run_Order = 1:nrow(temp.final))
				
				junk2 = data.frame(Factor = sfact1[temp.final$ID], Column = temp.final$Column)
				run.tab = table(junk2$Factor, junk2$Column)
			
			## assemble new sample names ##
			n.exp = gsub("[^[:alnum:]]","",exp.name)
			n.sfact1 = gsub("[^[:alnum:]]","",sfact1)
			n.stype = gsub("[^[:alnum:]]","",stype)
			n.sbio = gsub("[^[:alnum:]]","",sbio)

			if(omics.nm==TRUE){new.names = paste("OMICS",n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID], sprintf("%03d",temp.final$Overall_Prep_Order), sep="_")}else{
				new.names = paste(n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID],sprintf("%03d",temp.final$Overall_Prep_Order), sep="_")
			}

			## assemple results ##
			temp.final[,grep("ID", names(temp.final))] = new.names
			names(temp.final)[grep("ID", names(temp.final))] = "New_Name"
			names(temp.final)[sname.id] = "Original_Name"
				return(list(random = temp.final, tab1 = prep.tab, tab2 = run.tab))
		}
	}else{ if(rand.type == "a prep order only"){
			###### Prep Order #########
			# calculate batch sizes #
			numpfact.whole = table(sfact1) %/% num.prep.batches
			numpfact.part = table(sfact1) %% num.prep.batches

			samp.alloc = list()
			samp.unalloc = list()
			prep.batch.ind = list()
			for(i in 1:length(names(numpfact.whole))){
					temp = sampdata.wid[which(sfact1 == names(numpfact.whole[i])),]
				if(numpfact.whole[i] > 0 & numpfact.part[i] > 0){
					samp.id = sample(1:nrow(temp), num.prep.batches*numpfact.whole[i])
					prep.batch.ind[[i]] = rep(1:num.prep.batches, each = numpfact.whole[i])
					samp.alloc[[i]] = temp[samp.id,]
					samp.unalloc[[i]]= temp[-samp.id,]
					}
				if(numpfact.whole[i] == 0 & numpfact.part[i] > 0){
						prep.batch.ind[[i]] = NULL
						samp.alloc[[i]] = NULL
						samp.unalloc[[i]] = tem
					}
				if(numpfact.whole[i] > 0 & numpfact.part[i] == 0){
						samp.id = sample(1:nrow(temp), num.prep.batches*numpfact.whole[i])
						prep.batch.ind[[i]] = rep(1:num.prep.batches, each = numpfact.whole[i])
						samp.alloc[[i]] = temp[samp.id,]
						samp.unalloc[[i]] = NULL
						}
			}

			if(sum(numpfact.whole > 0) > 0 & sum(numpfact.part >0) > 0){
			samp.all = data.frame(do.call(rbind, samp.alloc), Prep_Batch = unlist(prep.batch.ind))
			samp.unall = cbind(do.call(rbind, samp.unalloc), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samp.unalloc))))
			}
			if(sum(numpfact.whole > 0) == 0 & sum(numpfact.part >0) > 0){
				samp.all = NULL
				samp.unall = cbind(do.call(rbind, samp.unalloc), Prep_Batch = rep(1:num.prep.batches, length = nrow(do.call(rbind, samp.unalloc))))
			}
			if(sum(numpfact.whole > 0) > 0 & sum(numpfact.part >0) == 0){
			samp.all = data.frame(do.call(rbind, samp.alloc), Prep_Batch = unlist(prep.batch.ind))
			samp.unall = NULL
			}

			samp.dat.wbatch = rbind(samp.all,samp.unall)

			temp.pord = list()
			for(j in 1:num.prep.batches){
				temp.dat = samp.dat.wbatch[which(samp.dat.wbatch$Prep_Batch==j),]
				temp.pord[[j]] = data.frame(temp.dat, Prep_Order = sample(1:nrow(temp.dat),nrow(temp.dat)))
			}

			samp.dat.wbatord = do.call(rbind, temp.pord)
			samp.dat.wbatord = samp.dat.wbatord[order(samp.dat.wbatord$Prep_Batch, samp.dat.wbatord$Prep_Order),]
			samp.dat.wbatord = data.frame(samp.dat.wbatord, Overall_Prep_Order = 1:nrow(samp.dat.wbatord))
	
				junk = data.frame(Factor = sfact1[samp.dat.wbatord$ID], Batch = samp.dat.wbatord$Prep_Batch)
				prep.tab = table(junk$Batch, junk$Factor)
	
				## assemble new sample names ##
				n.exp = gsub("[^[:alnum:]]","",exp.name)
				n.sfact1 = gsub("[^[:alnum:]]","",sfact1)
				n.stype = gsub("[^[:alnum:]]","",stype)
				n.sbio = gsub("[^[:alnum:]]","",sbio)

				if(omics.nm==TRUE){new.names = paste("OMICS",n.exp,n.sfact1[samp.dat.wbatord$ID], n.sbio[samp.dat.wbatord$ID], n.stype[samp.dat.wbatord$ID], sprintf("%03d",samp.dat.wbatch$Overall_Prep_Order), sep="_")}else{
					new.names = paste(n.exp,n.sfact1[samp.dat.wbatord$ID], n.sbio[samp.dat.wbatord$ID], n.stype[samp.dat.wbatord$ID],sprintf("%03d",samp.dat.wbatord$Overall_Prep_Order), sep="_")
				}

				## assemple results ##
				samp.dat.wbatord[,grep("ID", names(samp.dat.wbatord))] = new.names
				names(samp.dat.wbatord)[grep("ID", names(samp.dat.wbatord))] = "New_Name"
				names(samp.dat.wbatord)[sname.id] = "Old_Name"
				temp.final = samp.dat.wbatord
				
					return(list(random = temp.final, tab1 = prep.tab))
	}else { if(rand.type == "a run order only"){
			temp.d = sampdata.wid
			temp.fact = sfact1[temp.d$ID]

			numpf.whole = table(temp.fact) %/% num.columns
			numpf.part = table(temp.fact) %% num.columns

			temp.col.res.a = list()
			temp.col.res.u = list()
			run.col.ind = list()
			for(j in 1:length(names(numpf.whole))){
					tem = temp.d[which(temp.fact == names(numpf.whole[j])),]
					if(numpf.whole[j] > 0 & numpf.part[j] > 0){
					sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
					run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
					temp.col.res.a[[j]] = tem[sam.id,]
					temp.col.res.u[[j]]= tem[-sam.id,]
					}
					if(numpf.whole[j] > 0 & numpf.part[j] == 0){
					sam.id = sample(1:nrow(tem), num.columns*numpf.whole[j])
					run.col.ind[[j]] = rep(1:num.columns, each = numpf.whole[j])
					temp.col.res.a[[j]] = tem[sam.id,]
					temp.col.res.u[[j]]= NULL
					}
					if(numpf.whole[j] == 0 & numpf.part[j] > 0){
						run.col.ind[[j]] = NULL
						temp.col.res.a[[j]] = NULL
						temp.col.res.u[[j]] = tem
					}
			}

			if(sum(numpf.whole > 0) > 0 & sum(numpf.part > 0) > 0){
			temp.col.a = data.frame(do.call(rbind, temp.col.res.a), Column = unlist(run.col.ind))
			temp.col.u = data.frame(do.call(rbind, temp.col.res.u), Column = rep(1:num.columns, length = nrow(do.call(rbind, temp.col.res.u))))
			}
			if(sum(numpf.whole > 0) == 0 & sum(numpf.part > 0) > 0){
				temp.col.a = NULL
				temp.col.u = data.frame(do.call(rbind, temp.col.res.u), Column = rep(1:num.columns, length = nrow(do.call(rbind, temp.col.res.u))))
			}
			if(sum(numpf.whole > 0) > 0 & sum(numpf.part > 0) == 0){
				temp.col.a = data.frame(do.call(rbind, temp.col.res.a), Column = unlist(run.col.ind))
				temp.col.u = NULL
			}
		temp.run.res = rbind(temp.col.a, temp.col.u)

		temp.fina = list()
			for(j in 1:num.columns){
				tt = temp.run.res[which(temp.run.res$Column==j),]
				runord = sample(1:nrow(tt),nrow(tt))
				temp.fina[[j]] = data.frame(tt, Run_Order = runord)
				}

		temp.final = do.call(rbind, temp.fina)
		temp.final = temp.final[order(temp.final$Column,temp.final$Run_Order),]
		temp.final = data.frame(temp.final, Overall_Run_Order = 1:nrow(temp.final))

		junk2 = data.frame(Factor = sfact1[temp.final$ID], Column = temp.final$Column)
		run.tab = table(junk2$Factor, junk2$Column)

		## assemble new sample names ##
		n.exp = gsub("[^[:alnum:]]","",exp.name)
		n.sfact1 = gsub("[^[:alnum:]]","",sfact1)
		n.stype = gsub("[^[:alnum:]]","",stype)
		n.sbio = gsub("[^[:alnum:]]","",sbio)

		if(omics.nm==TRUE){new.names = paste("OMICS",n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID], sprintf("%03d",temp.final$Overall_Run_Order), sep="_")}else{
			new.names = paste(n.exp,n.sfact1[temp.final$ID], n.sbio[temp.final$ID], n.stype[temp.final$ID],sprintf("%03d",temp.final$Overall_Run_Order), sep="_")
		}

		## assemple results ##
		temp.final[,grep("ID", names(temp.final))] = new.names
		names(temp.final)[grep("ID", names(temp.final))] = "New_Name"
		names(temp.final)[sname.id] = "Old_Name"
		
			return(list(random = temp.final, tab1 = run.tab))
		
	}}
}

	
}