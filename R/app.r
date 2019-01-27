###############################################################################
# Globals
###############################################################################
#' is.ewas.object
#'
#' @param object a 'meffil' object (not an actual S3/4 R object)
is.ewas.object <- function(object) {
	is.list(object) && "class" %in% names(object) && object$class == "ewas"
}

#' scatter_thinning
#'
#' Remove points from a scatter plot where density is really high
#'
#' @param x x-coordinates vector
#' @param y y-coordinates vector
#' @param resolution number of partitions for the x and y-dimensions.
#' @param max_per_cell maximum number of points per x-y partition.
#' @return index into the points that omits points from x-y partitions
#' so that each has at most \code{max_per_cell} points.
scatter_thinning <- function(x, y, resolution = 100, max_per_cell = 100) {
	x.cell <- floor(
		(resolution - 1) * (x - min(x, na.rm = TRUE)) /
			diff(range(x, na.rm = TRUE))
	) + 1

	y.cell <- floor(
		(resolution - 1) * (y - min(y, na.rm = TRUE)) /
			diff(range(y, na.rm = TRUE))
	) + 1

	z.cell <- x.cell * resolution + y.cell
	frequency.table <- table(z.cell)
	frequency <- rep(0, max(z.cell, na.rm = TRUE))
	frequency[as.integer(names(frequency.table))] <- frequency.table
	f.cell <- frequency[z.cell]

	big.cells <- length(which(frequency > max_per_cell))

	sort(
		c(
			which(f.cell <= max_per_cell),
			sample(which(f.cell > max_per_cell),
						 size = big.cells * max_per_cell, replace = FALSE)
		),
		decreasing = FALSE
	)
}

#' manhattan_pre_proc_analysis

#' pre-procress and an analysis from an EWAS object in prepatation for a
#' manhattan plot
#' @param name
#' @param ewas
#' @param bi_dir

manhattan_pre_proc_analysis <- function(
	name, ewas, bi_dir = TRUE,
	chromosomes = paste("chr", c(1:22, "X", "Y"), sep="") # hard coding...
) {
	#expects: scatter_thinning
	chromosomes <- intersect(
		chromosomes,
		ewas$analyses[[name]]$table$chromosome
	)

	stats <- ewas$analyses[[name]]$table
	stats$chromosome <- factor(
		as.character(stats$chromosome),
		levels=chromosomes
	)
	stats$chr.colour <- 0
	# set every other chr colour to 1
	stats$chr.colour[
		stats$chromosome %in% chromosomes[seq(1, length(chromosomes), 2)]
		] <- 1
	p.values <- stats$p.value
	# set any p-values smaller than the smallest value of a double
	# allowed by the machine to the smallest value of a double allowed...~why?
	p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin

	stats$stat <- -log(p.values,10)
	if (bi_dir) {
		stats$stat <- stats$stat * sign(stats$coefficient)
	}
	stats <- stats[order(stats$stat, decreasing = TRUE),]

	chromosome.lengths <- sapply(
		chromosomes, function(chromosome) {
			max(stats$position[which(stats$chromosome == chromosome)])
		}
	)
	chromosome.lengths <- as.numeric(chromosome.lengths)

	chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
	names(chromosome.starts) <- c(chromosomes, "NA")
	# global/cumulative position
	stats$global <- stats$position + chromosome.starts[stats$chromosome] - 1

	selection.idx <- scatter_thinning(
		stats$global, stats$stat,
		resolution=100, max_per_cell=100
	)

	stats$probe <- rownames(stats)

	return(
		list(
			stats = stats[selection.idx,],
			name = name
		)
	)
}

#' manhattan_plotter
#'
#'
#' @param dat
#' @param threshold
#' @param prefix
#' @param plotly_prep
#' @param bi_dir
#'
manhattan_plotter <- function(
	dat,
	threshold = 1e-7,
	prefix = "Manhattan plot",
	plotly_prep = FALSE,
	bi_dir = TRUE,
	custom_colours = NULL#c("0"="paleblue","1"="darkblue","2"="red")
) {
	stats <- dat$stats
	name <- dat$name
	p <- ggplot(stats, aes(x = position, y = stat))

	if (bi_dir) {
		p <- p + geom_hline(yintercept = log(threshold,10), colour = "red")
	}

	p <- p +
		facet_grid(. ~ chromosome, space = "free_x", scales = "free_x") +
		theme(strip.text.x = element_text(angle = 90)) +
		labs(
			title = paste(prefix, ": ", name, sep = ""),
			x = "Position",
			y = bquote(-log[10]("p-value") * sign(beta))
		) +
		geom_hline(yintercept = -log(threshold,10), colour = "red") +
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

	if (!is.null(custom_colours)) {
		p <- p + scale_colour_manual(values = custom_colours) +
			guides(colour = TRUE)
	} else {
		p <- p + guides(colour = FALSE)
	}

	if (plotly_prep) {
		p <- p + geom_point(
			aes(
				colour = chr.colour,
				# tooltip contents
				text = paste0(
					"probe: ",probe,"\n",
					chromosome," : ",format(big.mark = ",",position),"\n",
					"-log10(p-value): ",sprintf("%.4g",abs(stat)),
					ifelse(
						"masks" %in% colnames(stats),
						paste0("\nmasks: ",masks),
						""
					)
				)
			)
		) +
			# ggplotly does not handle bquote well hence:
			labs(y = "-log10(p-value) x sign(beta)")
	} else {
		p <- p + geom_point(aes(colour = chr.colour))
	}
	return(p)
}

#' get_gene_cards
#'
#'
#' @param str
get_gene_cards <- function(str) {
	if (str != "") {
		genes <- str %>% strsplit(";") %>% unlist() %>% unique()
		lapply(genes,function(gene){
			#https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene
			paste0(
				"<a target='_blank' ",
				"href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
				gene,"'>",gene,"</a>"
			)
		}) %>% paste0(collapse = "; ")
	} else {"NA"}
}

###############################################################################
# UI
###############################################################################

#' Shiny app ui object
#'
#' @import shiny
#' @import shinydashboard
header <- dashboardHeader(title = "meffil EWAS viewer")

sidebar <- dashboardSidebar(
	sidebarMenu(
		menuItem("manhattan", tabName = "manhattan",icon = icon("chart-bar"))
	)
)

body <- dashboardBody(
	tabItems(
		tabItem(
			tabName = "manhattan",
			fluidRow(
				column(
					width = 4,
					fluidRow(
						uiOutput(outputId = "headlines")
					),
					fluidRow(
						box(
							width = 12,
							title = "Upload",
							collapsible = TRUE,
							status = "primary",
							solidHeader = TRUE,
							fileInput(
								inputId = "ewas",
								label="ewas object(s)",
								multiple = FALSE,
								accept = NULL, width = NULL,
								buttonLabel = "Browse...",
								placeholder = "No file selected"
							)
						)
					),
					fluidRow(
						box(
							title = "Inputs",
							status = "primary",
							width = 12,
							solidHeader = TRUE,
							collapsible = TRUE,
							column(
								width = 6,
								uiOutput(outputId = "pick_model"),
								numericInput(
									inputId = "alpha",
									label = "Significance Threshold",
									value = 1e-7,
									min = 0,
									max = 1
								),
								checkboxInput(
									inputId = "bi_dir",
									label = "bi-directional manhattan",
									value = FALSE
								)
							),
							column(
								width = 6,
								selectInput(
									inputId = "platform",
									label = "platform",
									choices = c("EPIC","450k"),
									selected = "EPIC"
								),
								#tableOutput(outputId = "selected_points"),
								numericInput(
									inputId = "DT_output_limit",
									label = "max rows to display",
									value = 500,
									min = 0,
									step = 1
								)
							)
						),
						box(
							width = 12,
							title = "Select A Point",
							solidHeader = TRUE,
							status = "info",
							uiOutput(outputId = "selected_points")
						)
					)
				),
				column(
					width = 8,
					box(
						title = "Manhattan",
						collapsible = TRUE,
						solidHeader = TRUE,
						status = "primary",
						width = 12,
						plotOutput(
							height = "60vh",
							outputId = "manhattan",
							click = "manhattan_click"#,
							# brush = brushOpts(
							# 	id = "manhattan_brush"
							# )

						)
					),
					box(
						title = "Probe Details",
						collapsible = TRUE,
						status = "primary",
						width = 12,
						solidHeader = TRUE,
						DT::dataTableOutput("manhattan_data")
					)
				)
			)
		)#
	)

)

ui <- dashboardPage(header, sidebar, body)

###############################################################################
# Server
###############################################################################

server <- function(input, output) {

	options(shiny.maxRequestSize=3000*1024^2)

	ewas <- reactive({
		if (is.null(input$ewas)){
			return(NULL)
		}
		readRDS(file = input$ewas$datapath)
	})

	preMan <- reactive({
		req(input$model, cancelOutput = TRUE)
		manhattan_pre_proc_analysis(
			input$model,
			bi_dir = input$bi_dir,
			ewas = ewas()
		)
	})

	observeEvent(input$ewas,{

		output$pick_model <- renderUI({
			selectInput(
				inputId = "model",
				label = "model",
				choices = names(ewas()$analyses),
				selected = names(ewas()$analyses)[1]
			)
		})

		output$headlines <- renderUI({
			req(input$model, cancelOutput = TRUE)
			box(
				title = "Outline",
				width = 12,
				color = "green",
				background = "green",
				collapsible = TRUE,
				solidHeader = TRUE,
				column(
					width = 6,
					h3(preMan()$stats %>%
						filter(p.value < input$alpha) %>%
						nrow()
					),
					renderText(
						paste0(
							"Significant sites (p < ",
							sprintf("%.6g",input$alpha),")"
						)
					),
					renderText(
						paste0(
							"Total sites: ",
							ewas()$analyses[[input$model]]$table %>%
								nrow() %>%
								format(big.mark = ",")
						)
					)
				),
				column(
					width = 6,
					h3("User Covariates"),
					renderText(
						ewas()$covariates %>%
							colnames() %>%
							paste0(collapse = " + ")
					)
				)
			)
		})
	})

	manifest <- reactive({
		switch(
			input$platform,
			"EPIC" = manifest_EPIC,
			"450k" = manifest_450k
		)
	})

	output$manhattan <- renderPlot({
		manhattan_plotter(
			preMan(),
			threshold = input$alpha,
			bi_dir = input$bi_dir
		)
	})

	# output$selected_points <- renderTable({
	# 	req(input$manhattan_brush,cancelOutput = TRUE)
	# 	# nearPoints(
	# 	# 	preMan()$stats,
	# 	# 	input$manhattan_click,
	# 	# 	addDist = FALSE,
	# 	# 	xvar = "position",
	# 	# 	yvar = "stat"
	# 	# ) %>%
	# 	brushedPoints(
	# 		preMan()$stats,
	# 		input$manhattan_brush,
	# 		xvar = "position",
	# 		yvar = "stat"
	# 	) %>%
	# 	select(chromosome,position,probe,p.value,fdr,coefficient)
	# })

	clickedPointRow <- reactive({
		req(input$manhattan_click,cancelOutput = TRUE)
		nearPoints(
			preMan()$stats,
			input$manhattan_click,
			addDist = FALSE,
			xvar = "position",
			yvar = "stat",
			maxpoints = 1
		)
	})

	#http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg19&position=chrX:Y-Z
	output$selected_points <- renderUI({
		pointDat <- clickedPointRow()
		expr = list(
			h4(paste0("Selected Probe: ",pointDat$probe)),
			HTML("<b>Probe Position:</b><br/>"),
			a(
				paste0(
					pointDat$chromosome," : ",
					format(pointDat$position,big.mark = ",")
				),
				href = paste0(
					"http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg19&position=",
					pointDat$chromosome,":",
					pointDat$position - 100,"-",
					pointDat$position +100
				),
				target="_blank"
			),
			HTML(
				paste0(
					" (hg19)<br/>",
					"<b>p-value:</b> ",sprintf("%.6g",pointDat$p.value),"<br/>",
					"<b>FDR:</b> ",sprintf("%.6g",pointDat$fdr),"<br/>",
					"<b>coefficient:</b> ",sprintf("%.6g",pointDat$coefficient),
					"<br/><b>Nearby Genes:</b> ",
					"<br/>"
				)
			),
			HTML(
				get_gene_cards(
					manifest() %>%
						filter(name == pointDat$probe) %>%
						pull(gene.symbol)
				)
			)
		)
	})

	output$manhattan_data <- DT::renderDataTable({

		# dplyr::select(chromosome,position,probe,p.value,fdr,coefficient) %>%
		# 	arrange(p.value) %>%
		# 	filter(p.value < input$alpha) %>%
		# 	head(n = input$DT_output_limit) %>%
		# 	DT::datatable()

		left_join(
			preMan()$stats %>%
				dplyr::select(probe,p.value,fdr,coefficient) %>%#chromosome,position,
				arrange(p.value) %>%
				filter(p.value < input$alpha) %>%
				head(n = input$DT_output_limit), #%>%
			manifest() %>%
				dplyr::rename(probe = name),
			by = "probe"
		) %>%
		dplyr::select(
			Chr=chromosome,Coord=position,probe,
			`p-value`=p.value,fdr,Coef=coefficient,
			`Gene Symbol`=gene.symbol#,masks
		) %>%
		DT::datatable(
			rownames = FALSE
		) %>%
		DT::formatSignif(c("p-value","fdr","Coef"), digits = 5) %>%
		DT::formatRound(c("Coord"),digits = 0)
	})

}
