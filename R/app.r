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

###############################################################################
# UI
###############################################################################

#' Shiny app ui object
#'
#' @import shiny

ui <- fluidPage(
	title = "EWAS results",
	sidebarLayout(
		sidebarPanel(
			h2("meffil EWAS viewer"),
			fileInput(
				inputId = "ewas", label="ewas object(s)", multiple = FALSE,
				accept = NULL, width = NULL,
				buttonLabel = "Browse...", placeholder = "No file selected"
			),
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
				value = TRUE
			),
			# h4("Selected Probes"),
			renderText("(Draw a box on the manhattan plot)"),
			uiOutput(outputId = "pick_model"),
			#tableOutput(outputId = "selected_points"),
			uiOutput(outputId = "selected_points"),
			numericInput(
				inputId = "DT_output_limit",
				label = "max rows to display",
				value = 500,
				min = 0,
				step = 1
			)
		),
		mainPanel(
			fluidRow(
				plotOutput(
					outputId = "manhattan",
					click = "manhattan_click"#,
					# brush = brushOpts(
					# 	id = "manhattan_brush"
					# )

				)
			),
			fluidRow(
				DT::dataTableOutput("manhattan_data")
			)
		)
	)
)

###############################################################################
# Server
###############################################################################

server <- function(input, output) {

	options(shiny.maxRequestSize=3000*1024^2)

	ewas <- reactive({
		#ewas <<- readRDS(file = input$ewas$datapath)
		if (is.null(input$ewas)){
			return(NULL)
		}
		readRDS(file = input$ewas$datapath)
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
	})

	preMan <- reactive({
		req(input$model, cancelOutput = TRUE)
		manhattan_pre_proc_analysis(
			input$model,
			bi_dir = input$bi_dir,
			ewas = ewas()
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
		#point <-
		nearPoints(
			preMan()$stats,
			input$manhattan_click,
			addDist = FALSE,
			xvar = "position",
			yvar = "stat",
			maxpoints = 1
		) #%>%
		#select(chromosome,position,probe,p.value,fdr,coefficient)
	})

	#http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg19&position=chrX:Y-Z
	output$selected_points <- renderUI({
		pointDat <- clickedPointRow()
		expr = list(
			h4(paste0("Selected Probe: ",pointDat$probe)),
			renderText("Probe Position: "),
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
				)
			),
			HTML(
				paste0(
					" (hg19)<br/>",
					"p-value: ",sprintf("%.6g",pointDat$p.value),"<br/>",
					"FDR: ",sprintf("%.6g",pointDat$fdr),"<br/>",
					"coefficient: ",sprintf("%.6g",pointDat$coefficient),
					"<br/><br/>"
				)
			)
		)
	})

	output$manhattan_data <- DT::renderDataTable({
		preMan()$stats %>%
			select(chromosome,position,probe,p.value,fdr,coefficient) %>%
			arrange(p.value) %>%
			filter(p.value < input$alpha) %>%
			head(n = input$DT_output_limit) %>%
			DT::datatable()
	})

}
