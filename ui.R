shinyUI(fluidPage(
	navbarPage("Omics Randomization Tool",
	
		tabPanel('Home',fluidPage(
			img(src = "bod.jpg", height = 300, width = 1500),
			h3("Omics Randomization Tool"),
			h4("Description"),
			p("The Omics Randomization Tool is setup to genrate randomized sample preparation and/or sample run order(s) based on a list of sample names uploaded by the user. The tool randomizes based on user inputs of experiemntal specifications (e.g. number of machine columns). The tool is currently capable of doing randomizations based on one and two factor experimental designs."),
			br(),
			h4("Requirements"),
			p("1. The uploaded list of sample names is in a file of .csv format."),
			p("2. One column of the .csv file should be dedicated to sample names. Other columns of information may be included and will be retained for the user."),
			p("3. Sample names contain relevant information (factor level, biological replicate, and sample type) in a character string with each component separated by the same delimiter (e.g. underscore, space, hyphen)"),
			p("4. The first row of the .csv file must be column names for any non-empty columns."),
			br(),
			h4("Instructions and Troubleshooting"),
			p("For help using the randomization tool click", strong(a("here", href = "http://prismwiki.pnl.gov/wiki/Omics_Randomization_Tool",target="_blank"))),
			p("For other questions, bug reports, etc., contact Lisa Bramer (lisa.bramer@pnnl.gov)")
			) # close fluid page
			
			), # close tabPanel
		
		## One Factor Experiment ##		
		tabPanel('One Factor Randomization', pageWithSidebar(
			# header #
			headerPanel(h3("One Factor Experiment")
			),
			
			# sidebar #
			sidebarPanel(
				fileInput('file1',h4('1. Upload List of Samples (.csv file)'), accept = c('text/csv','.csv')),
				uiOutput('colnames_select'),
				p(h4('Here is an example sample name:'), h5(textOutput('examp_name')	)),	
				p(h4("2. Specify Naming Structure")),
				h5("Based on the example, answer the following:"),
				uiOutput('experiment_selector_control'),
				conditionalPanel(condition = "input.experiment_selector == 'Other...'", textInput('man_expname','',value ="")),
				uiOutput('datatype_selector_control'),
				uiOutput('factor1_selector_control'),
				uiOutput('biorep_selector_control'),
				uiOutput('media_selector_control')
				),
			
			# main #
			mainPanel(
				tabsetPanel(id = 'subF1',
								tabPanel('Request',
								fluidPage(
									br(),
									fluidRow(
									column(6, 
										uiOutput('preprun_request_selector'),
										conditionalPanel(condition = "input.preprun_request == 'a prep order only' ",
										numericInput('numprep', 'How many prep batches are needed?', value = 1, min = 1, step = 1)),
										
										conditionalPanel(condition = "input.preprun_request == 'a run order only' ",
										numericInput('batchrun', 'How many run batches are needed?',value=1, min = 1, step = 1),
										numericInput('numrun', 'How many columns are on the instrument?', value = 1, min = 1, max = 4, step = 1)),
	
										
										conditionalPanel(condition = "input.preprun_request == 'prep and run orders' ",
										numericInput('numprep2', 'How many prep batches are needed?', value = 1, min = 1, step = 1),
										checkboxInput('preprun_connected','I need a run order that maintains prep batches'),
										conditionalPanel(condition = "input.preprun_connected == 0",
											numericInput('batchrun2', 'How many run batches are needed?', value = 1, min = 1, step = 1)),
										numericInput('numrun2', 'How many columns are on the instrument?', value = 1, min = 1, max = 4, step = 1)	
											)
										),
									column(6, 

										checkboxInput('paste_omics',"Paste an 'OMICS' prefix to my sample names"),
										actionButton("button","Run Randomization")
										)
										),
										h3("Randomization Results"),
							fluidRow(
								column(6,
								h4("Sample Prep Batch Distribution"),
								conditionalPanel(condition = "input.preprun_request != 'a run order only'", tableOutput('prep_tab_res'))	
								),
								column(6,
								h4("Sample Column Distribution"),
								conditionalPanel(condition = "input.preprun_request != 'a prep order only'",plotOutput('run_plot'))
									)	
									# close Row
									) 	
									
									) # close fluid page
									), # close tabPanel
								tabPanel('Randomization',
									br(),
									h3('Randomization Output'),
									br(),
									tableOutput('examp_res'),
									br(),
							downloadButton('downloadData', 'Download Randomization')
									) # close tabPanel
				)# close tabsetPanel
			
			
			) # close mainPanel
			) # close pagewithsidebar
), #close tabPanel
			## Two Factor Experiment ##		
tabPanel('Two Factor Randomization', pageWithSidebar(
	# header #
	headerPanel(h3("Two Factor Experiment")
	),
	
	# sidebar #
	sidebarPanel(
		fileInput('file2',h4('1. Upload List of Samples (.csv file)'), accept = c('text/csv','.csv')),
		uiOutput('colnames_select2'),
		p(h4('Here is an example sample name:'), h5(textOutput('examp_name2')	)),	
		p(h4("2. Specify Naming Structure")),
		h5("Based on the example, answer the following:"),
		uiOutput('experiment_selector_control2'),
		conditionalPanel(condition = "input.experiment_selector2 == 'Other...'", textInput('man_expname2','',value ="")),
		uiOutput('datatype_selector_control2'),
		uiOutput('factor1_selector_control2'),
		uiOutput('factor2_selector_control2'),
		uiOutput('biorep_selector_control2'),
		uiOutput('media_selector_control2')
		),
	
	# main #
	mainPanel(
			 tabsetPanel(id = 'subF2',
			 							tabPanel('Request',
			 							fluidPage(
			 								br(),
			 								fluidRow(
			 								column(6, 
			 									uiOutput('preprun_request_selector2'),
			 									conditionalPanel(condition = "input.preprun_request2 == 'a prep order only' ",
			 									numericInput('numprep2a', 'How many prep batches are needed?', value = 1, min = 1, step = 1)),
			 									
			 									conditionalPanel(condition = "input.preprun_request2 == 'a run order only' ",
			 									numericInput('batchrun2a', 'How many run batches are needed?',value=1, min = 1, step = 1),
			 									numericInput('numrun2a', 'How many columns are on the instrument?', value = 1, min = 1, max = 4, step = 1)),
			 
			 									
			 									conditionalPanel(condition = "input.preprun_request2 == 'prep and run orders' ",
			 									numericInput('numprep2b', 'How many prep batches are needed?', value = 1, min = 1, step = 1),
			 									checkboxInput('preprun_connected2','I need a run order that maintains prep batches'),
			 									conditionalPanel(condition = "input.preprun_connected2 == 0",
			 										numericInput('batchrun2b', 'How many run batches are needed?', value = 1, min = 1, step = 1)),
			 									numericInput('numrun2b', 'How many columns are on the instrument?', value = 1, min = 1, max = 4, step = 1)	
			 										)
			 									),
			 								column(6, 
			 								
			 									checkboxInput('paste_omics2',"Paste an 'OMICS' prefix to my sample names"),
			 									actionButton("button2","Run Randomization")
			 									)
			 									),
			 									h3("Randomization Results"),
			 						fluidRow(
			 							column(6,
			 							h4("Sample Prep Batch Distribution"),
			 							conditionalPanel(condition = "input.preprun_request != 'a run order only'", 
											h5("Factor1 Distribution by Batch"),
											tableOutput('prep_tab_res2a'),
											h5("Factor2 Distribution by Batch"),
											tableOutput('prep_tab_res2b')
											)	
			 							),
										column(6,
			 							h4("Sample Column Distribution"),
			 							conditionalPanel(condition = "input.preprun_request != 'a prep order only'",
											h5("Factor1 Distribution by Batch"),
											plotOutput('run_plot2a'),
											br(),
											h5("Factor2 Distribution by Batch"),
											plotOutput('run_plot2b'))
			 								)	
			 								# close Row
			 								) 	
			 								
											) # close fluid page
			 								), # close tabPanel
			 							tabPanel('Randomization',
			 								br(),
			 								h3('Randomization Output'),
			 								br(),
			 								tableOutput('examp_res2'),
			 								br(),
			 						downloadButton('downloadData2', 'Download Randomization')
			 								) # close tabPanel
			 			)# close tabsetPanel
			 		)
		)#close pagewithsidebar
		)# close tabPanel
)#close navbar page
)#close fluid page
)#close shinyUI			