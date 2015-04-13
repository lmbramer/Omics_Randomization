library(shiny)
library(ggplot2)
source("find_delimiter.R", chdir = TRUE)
source("onefact.R", chdir = TRUE)
source("twofact.R", chdir = TRUE)
source("twofactw.R", chdir = TRUE)

shinyServer(function(input, output){

###################################### One Factor ##################################
	
	## read in data from csv ##
	## then output data and column names ##
	dataInput <- reactive({
		if (is.null(input$file1)){return(list(dat = NULL, nm.choices = NULL))}else{
				data = read.csv(input$file1$datapath, header=T)
				nm.choices = names(data)
				return(list(dat = data, nam.choices = nm.choices))
			}
		})
		
	
########### Data Formatting ################
	
	## have user select column for sample names ##	
		colnames_selector_control = reactive({
			if(is.null(dataInput()$dat)){nm_choices = c("")}else{nm_choices = dataInput()$nam.choices}
			nm_choices
		})
		
		## UI output for selecting column with sample names ##
		output$colnames_select = renderUI({
		selectInput('colnames_selector','Select the name of the column that contains sample names:', colnames_selector_control(), selectize=F)
			
		})
		
			
	## 	pull example name for display ##
		output$examp_name = renderText({
			if(is.null(dataInput()$dat) ){""}else{
			as.character(dataInput()$dat[1,which( dataInput()$nam.choices == input$colnames_selector )])
			}
		})


		delim_select = reactive({
			if(input$colnames_selector==""){""}else{
				find_delimiter(as.character(dataInput()$dat[1,input$colnames_selector]))
			}
		})
			
			
	## generate list of string parts ##	
		selector_variables = reactive({
			if(is.null(dataInput()$dat)){""}else{
					nm.col = which(dataInput()$nam.choices == input$colnames_selector)
					unlist(strsplit(as.character(dataInput()$dat[1,nm.col]),delim_select()))
				}
				
			})


	## UI output for experiment selector ##
		output$experiment_selector_control = renderUI({
			selectInput('experiment_selector','Select the field that indicates the experiment name:', c(selector_variables(),"Other..."), selectize=F)	
		})
		

	## UI output for sample type selector ##
			output$datatype_selector_control = renderUI({
			selectInput('datatype_selector','Select the field that indicates the sample type:',	selector_variables(), selectize=F)
			})


	## UI output for factor 1 selector ##	
			output$factor1_selector_control = renderUI({
				selectInput('factor1_selector','Select the field that indicates your first factor:',selector_variables(), selectize=F)
			})


	## UI output for biorep selector ##	
			output$biorep_selector_control = renderUI({
				selectInput('biorep_selector','Select the field that indicates the biological replicate:',selector_variables(), selectize=F)
			})



	## pull sample types that exist in the data ##
			media_types = reactive({
				if(is.null(dataInput()$dat)){""}else{
					col.id = which( dataInput()$nam.choices == input$colnames_selector )
					me.id = which(selector_variables() == input$datatype_selector)
					unique(unlist(lapply(strsplit(as.character(dataInput()$dat[,col.id]), delim_select()), function(x) x[me.id])))
				}
			})

	## UI output for media type selector ##
	output$media_selector_control = renderUI({
				selectInput('media_selector','Select sample type to randomize',media_types(),selectize=F)

			})



############ Prep and Run Randomization Details ############

	output$preprun_request_selector = renderUI({
		selectInput('preprun_request',h4('I need randomizations for . . .'), c("prep and run orders", "a prep order only","a run order only"),selectize=F)
	})
	
	num.ids = reactive({
		if(is.null(dataInput()$dat)){
		list(col.id = NULL, me.id = NULL, bio.id = NULL, bio.id = NULL, fact.id = NULL, exp.n = NULL)
		}else{
				col.id = which( dataInput()$nam.choices == input$colnames_selector )
				me.id = which(selector_variables() == input$datatype_selector)
				bio.id = which(selector_variables() == input$biorep_selector)
				fact.id = which(selector_variables() == input$factor1_selector)
				if(input$experiment_selector== 'Other...'){exp.n = input$man_expname}else{
					exp.n = input$experiment_selector
				}
				list(col.id = col.id, me.id = me.id, bio.id = bio.id, fact.id = fact.id, exp.n = exp.n)
				}	
		})

########## Perform Randomization ############
	preprun_results_go = eventReactive(input$button, ignoreNULL=TRUE,{
	
		if(input$preprun_request=="prep and run orders" & input$preprun_connected==TRUE){
			res = onefact_prep_run_connected(samp.data = dataInput()$dat,
								name.delim = delim_select(),
								num.prep.batches = input$numprep2,
								num.columns = input$numrun2,
								exp.name = num.ids()$exp.n,
								sname.id = num.ids()$col.id,
								stype.id = num.ids()$me.id,
								sbio.id = num.ids()$bio.id,
								sfact1.id = num.ids()$fact.id,
								samptype.rand = input$media_selector,
								omics.nm = input$paste_omics)
								
						rets =	list(random = res$rand.data, tabrun = res$run.tab, tabprep = res$prep.tab)
		}
		if(input$preprun_request=="prep and run orders" & input$preprun_connected==FALSE){
				res = onefact_prep_run_notconnected(samp.data = dataInput()$dat,
									name.delim = delim_select(),
									num.prep.batches = input$numprep2,
									num.run.batches = input$batchrun2,
									num.columns = input$numrun2,
									exp.name = num.ids()$exp.n,
									sname.id = num.ids()$col.id,
									stype.id = num.ids()$me.id,
									sbio.id = num.ids()$bio.id,
									sfact1.id = num.ids()$fact.id,
									samptype.rand = input$media_selector,
									omics.nm = input$paste_omics)
							rets =	list(random = res$rand.data, tabrun = res$run.tab, tabprep = res$prep.tab)}
		if(input$preprun_request=="a prep order only"){
					res = onefact_prep(samp.data = dataInput()$dat,
										name.delim = delim_select(),
										num.prep.batches = input$numprep,
										exp.name = num.ids()$exp.n,
										sname.id = num.ids()$col.id,
										stype.id = num.ids()$me.id,
										sbio.id = num.ids()$bio.id,
										sfact1.id = num.ids()$fact.id,
										samptype.rand = input$media_selector,
										omics.nm = input$paste_omics)
								rets =	list(random = res$rand.data, tabprep = res$prep.tab, tabrun = NULL)}
		if(input$preprun_request=="a run order only"){
					res = onefact_run(samp.data = dataInput()$dat,
										name.delim = delim_select(),
										num.run.batches = input$batchrun,
										num.columns = input$numrun,
										exp.name = num.ids()$exp.n,
										sname.id = num.ids()$col.id,
										stype.id = num.ids()$me.id,
										sbio.id = num.ids()$bio.id,
										sfact1.id = num.ids()$fact.id,
										samptype.rand = input$media_selector,
										omics.nm = input$paste_omics)
								 rets = list(random = res$rand.data, tabrun = res$run.tab, tabprep = NULL)
								}
								rets
							
	}) 
	
################ Output Results ##################
	
out_cond = reactive({ 

	if(is.null(dataInput()$dat)){NULL}else{
	if(input$preprun_request=="prep and run orders" & input$preprun_connected==FALSE){
		"cond1"}
	if(input$preprun_request=="prep and run orders" & input$preprun_connected==TRUE){
			"cond2"}	
	if(input$preprun_request=="a prep order only" & input$preprun_connected==TRUE){
			"cond3"}
	if(input$preprun_request=="a run order only" & input$preprun_connected==TRUE){
			"cond4"}	
	}
	})



output$examp_res = renderTable({
	head(preprun_results_go()$random)
	},include.rownames=FALSE)
	

output$prep_tab_res = renderTable({
	
		preprun_results_go()$tabprep
		})
		
		
output$run_plot = renderPlot({

			if(max(as.numeric(preprun_results_go()$tabrun[,1])) > 1){
				nm.bat = max(as.numeric(preprun_results_go()$tabrun[,1]))
				dat = preprun_results_go()$tabrun
				levels(dat$Run_Batch) = paste("Batch",1:nm.bat, sep = " ")
					ggplot(data = dat, aes(x = factor(Factor1), y = Freq, fill = factor(Block))) +
					geom_bar(stat="identity",position="dodge") +
					facet_wrap(~Run_Batch) +
					xlab("Factor") +
					ylab("Number of Samples") +
					theme_bw() +
					scale_fill_discrete(guide = guide_legend(title = "Column/Block")) }else{
						ggplot(data = preprun_results_go()$tabrun, aes(x = factor(Factor1), y = Freq, fill = factor(Block))) +
						geom_bar(stat="identity",position="dodge") +
						xlab("Factor") +
						ylab("Number of Samples") +
						theme_bw() +
						scale_fill_discrete(guide = guide_legend(title = "Column/Block"))
					}
				

	})
	
	output$downloadData <- downloadHandler(
	   filename = function() { paste(num.ids()$exp.n,'randomization', '.csv', sep='') },
	    content = function(filename) {
	    write.csv(preprun_results_go()$random, filename,row.names=FALSE)
	    }
	  )


	###################################### Two Factor ##################################

	 	## read in data from csv ##
	 	## then output data and column names ##
	 	dataInput2 <- reactive({
	 		if (is.null(input$file2)){return(list(dat = NULL, nm.choices = NULL))}else{
	 				data2 = read.csv(input$file2$datapath, header=T)
					nm.choices2 = names(data2)
	 				return(list(dat = data2, nam.choices = nm.choices2))
				}
	 		})
	 
	 
	 ########### Data Formatting ################
	 
	 	## have user select column for sample names ##	
			colnames_selector_control2 = reactive({
	 			if(is.null(dataInput2()$dat)){nm_choices2 = c("")}else{nm_choices2 = dataInput2()$nam.choices}
	 			nm_choices2
	 		})
	 
	 		## UI output for selecting column with sample names ##
	 		output$colnames_select2 = renderUI({
	 		selectInput('colnames_selector2','Select the name of the column that contains sample names:', colnames_selector_control2(), selectize=F)
	
	 		})
	 
	 
	 	## 	pull example name for display ##
			output$examp_name2 = renderText({
				if(is.null(dataInput2()$dat) ){""}else{
				as.character(dataInput2()$dat[1,which( dataInput2()$nam.choices == input$colnames_selector2 )])
				}
			})
	
	 
	 		delim_select2 = reactive({
	 			if(input$colnames_selector2==""){""}else{
	 				find_delimiter(as.character(dataInput2()$dat[1,input$colnames_selector2]))
	 			}
	 		})
	 
	 
	 	## generate list of string parts ##	
	 		selector_variables2 = reactive({
	 			if(is.null(dataInput2()$dat)){""}else{
	 					nm.col = which(dataInput2()$nam.choices == input$colnames_selector2)
	 					unlist(strsplit(as.character(dataInput2()$dat[1,nm.col]),delim_select2()))
	 				}
	 
	 			})
	 
	 
	 	## UI output for experiment selector ##
	 		output$experiment_selector_control2 = renderUI({
	 			selectInput('experiment_selector2','Select the field that indicates the experiment name:', c(selector_variables2(),"Other..."), selectize=F)	
	 		})
	 
	 
	 	## UI output for sample type selector ##
				output$datatype_selector_control2 = renderUI({
	 			selectInput('datatype_selector2','Select the field that indicates the sample type:',	selector_variables2(), selectize=F)
	 			})
	 
	 
	 	## UI output for factor 1 selector ##	
	 			output$factor1_selector_control2 = renderUI({
	 				selectInput('factor1_selector2','Select the field that indicates your first factor:',selector_variables2(), selectize=F)
	 			})
	 

		## UI output for factor 1 selector ##	
			 	output$factor2_selector_control2 = renderUI({
			 			selectInput('factor2_selector2','Select the field that indicates your second factor:',selector_variables2(), selectize=F)
			 			})
	 
		## UI output for biorep selector ##	
	 			output$biorep_selector_control2 = renderUI({
	 				selectInput('biorep_selector2','Select the field that indicates the biological replicate:',selector_variables2(), selectize=F)
	 			})
	 
	 
	 
	 	## pull sample types that exist in the data ##
	 			media_types2 = reactive({
	 				if(is.null(dataInput2()$dat)){""}else{
	 					col.id = which( dataInput2()$nam.choices == input$colnames_selector2 )
	 					me.id = which(selector_variables2() == input$datatype_selector2)
	 					unique(unlist(lapply(strsplit(as.character(dataInput2()$dat[,col.id]), delim_select2()), function(x) x[me.id])))
	 				}
	 			})
	 
	 	## UI output for media type selector ##
	 	output$media_selector_control2 = renderUI({
	 				selectInput('media_selector2','Select sample type to randomize',media_types2(),selectize=F)
	
	 			})
	 
	 
	 
	 ############ Prep and Run Randomization Details ############
	 
	 	output$preprun_request_selector2 = renderUI({
	 		selectInput('preprun_request2',h4('I need randomizations for . . .'), c("prep and run orders", "a prep order only","a run order only"),selectize=F)
	 	})
	 
	 	num.ids2 = reactive({
	 		if(is.null(dataInput2()$dat)){
	 		list(col.id = NULL, me.id = NULL, bio.id = NULL, bio.id = NULL, fact.id = NULL, exp.n = NULL)
	 		}else{
	 				col.id = which( dataInput2()$nam.choices == input$colnames_selector2 )
	 				me.id = which(selector_variables2() == input$datatype_selector2)
	 				bio.id = which(selector_variables2() == input$biorep_selector2)
	 				fact1.id = which(selector_variables2() == input$factor1_selector2)
					fact2.id = which(selector_variables2() == input$factor2_selector2)
	 				if(input$experiment_selector2== 'Other...'){exp.n = input$man_expname2}else{
	 					exp.n = input$experiment_selector2
	 				}
	 				list(col.id = col.id, me.id = me.id, bio.id = bio.id, fact1.id = fact1.id, fact2.id = fact2.id,exp.n = exp.n)
	 				}	
	 		})
	 
	 ########## Perform Randomization ############
	 	preprun_results_go2 = eventReactive(input$button2, ignoreNULL=TRUE,{
	
	 		if(input$preprun_request2=="prep and run orders" & input$preprun_connected2==TRUE){
				res = twofact_prep_run_connectedw(samp.data = dataInput2()$dat,
									name.delim = delim_select2(),
	 								num.prep.batches = input$numprep2b,
									num.columns = input$numrun2b,
	 								exp.name = num.ids2()$exp.n,
	 								sname.id = num.ids2()$col.id,
	 								stype.id = num.ids2()$me.id,
	 								sbio.id = num.ids2()$bio.id,
	 								sfact1.id = num.ids2()$fact1.id,
									sfact2.id = num.ids2()$fact2.id,
	 								samptype.rand = input$media_selector2,
	 								omics.nm = input$paste_omics2)
	 
	 				rets =	list(random = res$rand.data, tabrun = res$run.tab, tabprep1 = res$prep.tab1, tabprep2 = res$prep.tab2)
	 		}
	 		if(input$preprun_request2=="prep and run orders" & input$preprun_connected2==FALSE){
	 				res = twofact_prep_run_notconnectedw(samp.data = dataInput2()$dat,
	 									name.delim = delim_select2(),
	 									num.prep.batches = input$numprep2b,
	 									num.run.batches = input$batchrun2b,
	 									num.columns = input$numrun2b,
	 									exp.name = num.ids2()$exp.n,
	 									sname.id = num.ids2()$col.id,
	 									stype.id = num.ids2()$me.id,
	 									sbio.id = num.ids2()$bio.id,
	 									sfact1.id = num.ids2()$fact1.id,
										sfact2.id = num.ids2()$fact2.id,
	 									samptype.rand = input$media_selector2,
	 									omics.nm = input$paste_omics2)
	 				rets =	list(random = res$rand.data, tabrun = res$run.tab, tabprep1 = res$prep.tab1, tabprep2 = res$prep.tab2)}
	 		if(input$preprun_request2=="a prep order only"){
	 					res = twofact_prepw(samp.data = dataInput2()$dat,
	 										name.delim = delim_select2(),
	 										num.prep.batches = input$numprep2a,
	 										exp.name = num.ids2()$exp.n,
	 										sname.id = num.ids2()$col.id,
	 										stype.id = num.ids2()$me.id,
	 										sbio.id = num.ids2()$bio.id,
	 										sfact1.id = num.ids2()$fact1.id,
											sfact2.id = num.ids2()$fact2.id,
	 										samptype.rand = input$media_selector2,
	 										omics.nm = input$paste_omics2)
	 				rets =	list(random = res$rand.data, tabprep1 = res$prep.tab1, tabprep2 = res$prep.tab2, tabtabrun = NULL)}
			if(input$preprun_request2=="a run order only"){
	 					res = twofact_run(samp.data = dataInput2()$dat,
											name.delim = delim_select2(),
	 										num.run.batches = input$batchrun2a,
	 										num.columns = input$numrun2a,
	 										exp.name = num.ids2()$exp.n,
	 										sname.id = num.ids2()$col.id,
	 										stype.id = num.ids2()$me.id,
	 										sbio.id = num.ids2()$bio.id,
	 										sfact1.id = num.ids2()$fact1.id,
											sfact2.id = num.ids2()$fact2.id,
	 										samptype.rand = input$media_selector2,
	 										omics.nm = input$paste_omics2)
	 								 rets = list(random = res$rand.data, tabrun = res$run.tab, tabprep1 = NULL, tabprep2 = NULL)
	 								}
	 								rets
	 
		}) 
	
	 ################ Output Results ##################
	 
	 out_cond2 = reactive({ 
	 
	 	if(is.null(dataInput2()$dat)){NULL}else{
	 	if(input$preprun_request2=="prep and run orders" & input$preprun_connected2==FALSE){
	 		"cond1"}
	 	if(input$preprun_request2=="prep and run orders" & input$preprun_connected2==TRUE){
	 			"cond2"}	
	 	if(input$preprun_request2=="a prep order only" & input$preprun_connected2==TRUE){
	 			"cond3"}
	 	if(input$preprun_request2=="a run order only" & input$preprun_connected2==TRUE){
	 			"cond4"}	
	 	}
	 	})
	 
	 
	 
	 output$examp_res2 = renderTable({
	 	head(preprun_results_go2()$random)
	 	},include.rownames=FALSE)
	 
	 
	 output$prep_tab_res2a = renderTable({
	  		preprun_results_go2()$tabprep1
	 		})
	 
	output$prep_tab_res2b = renderTable({
			 preprun_results_go2()$tabprep2
			 })
	 
	 output$run_plot2a = renderPlot({
	 
	 			if(max(as.numeric(preprun_results_go2()$tabrun[,1])) > 1){
	 				nm.bat = max(as.numeric(preprun_results_go2()$tabrun[,1]))
	 				dat = preprun_results_go2()$tabrun
	 				levels(dat$Run_Batch) = paste("Batch",1:nm.bat, sep = " ")
	 					g1 = ggplot(data = dat, aes(x = factor(Factor1), y = Freq, fill = factor(Block))) +
	 					geom_bar(stat="identity",position="dodge") +
	 					facet_wrap(~Run_Batch, ncol=2) +
	 					xlab("Factor") +
	 					ylab("Number of Samples") +
	 					theme_bw() +
	 					scale_fill_discrete(guide = guide_legend(title = "Column/Block")) 
	
					
	}else{
	 					g1 = ggplot(data = preprun_results_go2()$tabrun, aes(x = factor(Factor1), y = Freq, fill = factor(Block))) +
	 						geom_bar(stat="identity",position="dodge") +
	 						xlab("Factor") +
	 						ylab("Number of Samples") +
	 						theme_bw() +
	 						scale_fill_discrete(guide = guide_legend(title = "Column/Block"))
	
	 					}
	 
	 		g1
	 	})
	 
		 output$run_plot2b = renderPlot({

		 			if(max(as.numeric(preprun_results_go2()$tabrun[,1])) > 1){
		 				nm.bat = max(as.numeric(preprun_results_go2()$tabrun[,1]))
		 				dat = preprun_results_go2()$tabrun
		 				levels(dat$Run_Batch) = paste("Batch",1:nm.bat, sep = " ")
		 					g1 = ggplot(data = dat, aes(x = factor(Factor2), y = Freq, fill = factor(Block))) +
		 					geom_bar(stat="identity",position="dodge") +
		 					facet_wrap(~Run_Batch, ncol=2) +
		 					xlab("Factor") +
		 					ylab("Number of Samples") +
		 					theme_bw() +
		 					scale_fill_discrete(guide = guide_legend(title = "Column/Block")) 


		}else{
		 					g1 = ggplot(data = preprun_results_go2()$tabrun, aes(x = factor(Factor2), y = Freq, fill = factor(Block))) +
		 						geom_bar(stat="identity",position="dodge") +
		 						xlab("Factor") +
		 						ylab("Number of Samples") +
		 						theme_bw() +
		 						scale_fill_discrete(guide = guide_legend(title = "Column/Block"))

		 					}

		 		g1
		 	})
		
	 	output$downloadData2 <- downloadHandler(
	 	   filename = function() { paste(num.ids2()$exp.n,'randomization', '.csv', sep='') },
	 	    content = function(filename) {
	 	    write.csv(preprun_results_go2()$random, filename,row.names=FALSE)
	 	    }
	 	  )
	 
	 

}
)