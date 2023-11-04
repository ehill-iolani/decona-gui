library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(stringr)
library(dplyr)

# Define accessory fucnction
# Checks if an integer is of 0 length
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

# UI setup
ui <- shinyUI(fluidPage(
  titlePanel("AIN Decona_plus"),
  sidebarLayout(
    sidebarPanel(
      helpText("Select a directory where the barcode directories containing the fastq or fastq.gz files are located"),
      fileInput("fastqs", "Input fastq file directories", accept = c("/"), multiple = TRUE),
      fileInput("database", "Upload database .fasta file", accept = c("fasta")),
      numericInput("lowerlength", "Minimum amplicon length", value = 170),
      numericInput("upperlength", "Maximum amplicon length", value = 230),
      numericInput("quality", "Minimum Quality Score", value = 10),
      numericInput("clusterid", "Cluster Percent Identity", value = 0.95),
      numericInput("minclustersize", "Minimum Cluster Size", value = 5),
      numericInput("kmer", "Kmer-length", value = 10),
      numericInput("threads", "Threads", value = 2),
      actionButton("run", "Run Decona!")),
    mainPanel(
      verbatimTextOutput("results"))
  )
))

# Server setup
server <- function(input, output) {
  # When the run button is clicked, execute decona
  observeEvent(input$run, {
    showModal(modalDialog("Running the Decona pipeline, please be patient...", footer = NULL))

    # Set large upload size limit (server side)
    options(shiny.maxRequestSize = 70 * 1024^2)

    # Make a directory to store input files
    dir.create(file.path("/home/processing"))

    # Make a directory for each barcode
    files <- input$fastqs$name
    bc_names <- str_extract_all(input$fastqs$name, "barcode[0-9]+") %>% unique()
    bc_names <- unlist(bc_names)
    for (i in seq_along(bc_names)) {
      dir.create(file.path("/home/processing", bc_names[i]))
    }

    # Write fastq files to their respective barcode directories
    for (i in seq_along(bc_names)) {
      bucket <- which(str_detect(files, bc_names[i]))
      file.copy(input$fastqs$datapath[bucket], file.path("/home/processing", bc_names[i]))
    }

    # Move the database file to the input directory
    file.copy(input$database$datapath, file.path("/home/processing"))

    # Set the temporary working directory
    setwd("/home/processing")

    # prep the decona command
    decona_command <- paste0("conda run -n decona decona -f -l ", input$lowerlength, " -m ", input$upperlength, " -q ", input$quality, " -c ", input$clusterid, " -n ", input$minclustersize, " -k ", input$kmer, " -T ", input$threads, " -B 0.fasta")

    # combine the docker and decona commands
    command <- paste(decona_command)

    # run the decona pipeline
    system(command, intern = TRUE)

    # Display the BLAST results of the cluters
    output$results <- renderPrint({
      # Read in the raw text output
      dat <- read.delim(file.path("/home/processing/result/Racon", "summary_BLAST_out_racon_clusters_barcode01_concatenated.txt"), sep = "\t", header = FALSE, skip = 1, col.names = paste0("V",seq_len(100)), fill = TRUE)
      check <- grep("TRUE", is.na(dat$V1))

      # If there are no clusters, return an empty data frame
      if(is.integer0(check) == TRUE) {
        sdat <- dat[, c(1, 4, 8, 11, 12)]
      } else {
        dat <- dat[-c(grep("TRUE", is.na(dat$V1))), ]
        sdat <- dat[, c(1, 4, 8, 11, 12)]}

      # Rename the columns
      names(sdat) <- c("reads", "percent_id", "e_value", "genus", "species")
      sdat$percent_id[grep("TRUE", is.na(sdat$percent_id))] <- 0
      sdat
    })

    removeModal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
