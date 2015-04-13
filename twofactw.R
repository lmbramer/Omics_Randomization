twofact_prep_run_connectedw <- function(samp.data, name.delim, num.prep.batches, num.columns, exp.name, sname.id, stype.id, sbio.id, sfact1.id, sfact2.id, samptype.rand, omics.nm){

alt1 = twofact_prep_run_connected(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, num.columns = num.columns, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact1.id, sfact2.id = sfact2.id, samptype.rand = samptype.rand, omics.nm = omics.nm)

alt2 = twofact_prep_run_connected(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, num.columns = num.columns, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact2.id, sfact2.id = sfact1.id, samptype.rand = samptype.rand, omics.nm = omics.nm)


clac <- function(x){
	cnt = 1
	temp = matrix(0,nrow(x),(ncol(x))^2)
	for(i in 1:ncol(x)){
		for(j in 1:ncol(x)){
			temp[,cnt] = x[,i] - x[,j]
		cnt = cnt + 1
		}
	}
max(apply(temp, 1, max))	
}

a1 = max(clac(alt1$prep.tab1),clac(alt1$prep.tab2))
a2 = max(clac(alt2$prep.tab1),clac(alt2$prep.tab2))

if(a1 <= a2){return(alt1)}else{return(alt2)}

}


twofact_prep_run_notconnectedw <- function(samp.data, name.delim, num.prep.batches, num.run.batches, num.columns, exp.name, sname.id, stype.id, sbio.id, sfact1.id, sfact2.id, samptype.rand, omics.nm){

alt1 = twofact_prep_run_notconnected(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, num.run.batches = num.run.batches, num.columns = num.columns, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact1.id, sfact2.id = sfact2.id, samptype.rand = samptype.rand, omics.nm = omics.nm)

alt2 = twofact_prep_run_notconnected(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, num.run.batches = num.run.batches, num.columns = num.columns, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact2.id, sfact2.id = sfact1.id, samptype.rand = samptype.rand, omics.nm = omics.nm)


clac <- function(x){
	cnt = 1
	temp = matrix(0,nrow(x),(ncol(x))^2)
	for(i in 1:ncol(x)){
		for(j in 1:ncol(x)){
			temp[,cnt] = x[,i] - x[,j]
		cnt = cnt + 1
		}
	}
max(apply(temp, 1, max))	
}

a1 = max(clac(alt1$prep.tab1),clac(alt1$prep.tab2))
a2 = max(clac(alt2$prep.tab1),clac(alt2$prep.tab2))

if(a1 <= a2){return(alt1)}else{return(alt2)}

}

twofact_prepw <- function(samp.data, name.delim, num.prep.batches, exp.name, sname.id, stype.id, sbio.id, sfact1.id, sfact2.id, samptype.rand, omics.nm){

alt1 = twofact_prep(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact1.id, sfact2.id = sfact2.id, samptype.rand = samptype.rand, omics.nm = omics.nm)

alt2 = twofact_prep(samp.data = samp.data, name.delim = name.delim, num.prep.batches = num.prep.batches, exp.name = exp.name, sname.id = sname.id, stype.id = stype.id, sbio.id = sbio.id, sfact1.id = sfact2.id, sfact2.id = sfact1.id, samptype.rand = samptype.rand, omics.nm = omics.nm)


clac <- function(x){
	cnt = 1
	temp = matrix(0,nrow(x),(ncol(x))^2)
	for(i in 1:ncol(x)){
		for(j in 1:ncol(x)){
			temp[,cnt] = x[,i] - x[,j]
		cnt = cnt + 1
		}
	}
max(apply(temp, 1, max))	
}

a1 = max(clac(alt1$prep.tab1),clac(alt1$prep.tab2))
a2 = max(clac(alt2$prep.tab1),clac(alt2$prep.tab2))

if(a1 <= a2){return(alt1)}else{return(alt2)}

}