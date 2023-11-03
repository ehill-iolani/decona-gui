library(shiny)
library(shinyFiles)
library(shinycssloaders)

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
      shinyDirButton("dir", "Input directory", "Upload"),
      verbatimTextOutput("dir", placeholder = TRUE),
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
  # Parses the input directory
  shinyDirChoose(
    input,
    "dir",
    roots = c(home = "~")
    #debatable if we want to allow the user to select the file type....
    #i think im just gonna assume theyre not braindead and can
    #read the documenaion to see what file types are allowed
    #filetypes = c("fastq", "fastq.gz")
  )

  # Creates a reactive value for the input directory
  global <- reactiveValues(datapath = getwd())
  dir <- reactive(input$dir)

  # Outputs the input directory path
  output$dir <- renderText({
    global$datapath
  })

  # Updates the input directory path
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })

  # When the run button is clicked, execute decona
  observeEvent(input$run, {
    showModal(modalDialog("Running the Decona pipeline, please be patient...", footer = NULL))
    files <- list.files(global$datapath)

    # Make a directory to store input files
    dir.create(file.path(global$datapath, "processing"))
    # Move the database file to the input directory
    file.copy(input$database$datapath, file.path(global$datapath, "processing"))
    # Move the fastq files to the input directory
    file.copy(file.path(global$datapath, files), file.path(global$datapath, "processing"), recursive = TRUE)

    ### IMPORTANT: THIS WILL NEED TO BE CHANGED WHEN DEPLOYED B/C IT WILL ALL BE CONTAINED IN THE DOCKER IMAGE
    # prep the docker command
    docker_command <- paste0("docker run --name=decona_app ", "-v ", global$datapath, "/processing/", ":/home/data ", "decona:dev")

    # prep the decona command
    decona_command <- paste0("decona -f -l ", input$lowerlength, " -m ", input$upperlength, " -q ", input$quality, " -c ", input$clusterid, " -n ", input$minclustersize, " -k ", input$kmer, " -T ", input$threads, " -B 0.fasta")

    # combine the docker and decona commands
    command <- paste(docker_command, decona_command)

    # run the decona pipeline
    system(command, intern = TRUE)

    # delete the docker container
    system("docker rm decona_app", intern = TRUE)

    # Display the BLAST results of the cluters
    output$results <- renderPrint({
      # Read in the raw text output
      dat <- read.delim(file.path(global$datapath, "processing/result/Racon", "summary_BLAST_out_racon_clusters_barcode01_concatenated.txt"), sep = "\t", header = FALSE, skip = 1, col.names = paste0("V",seq_len(100)), fill = TRUE)
      check <- grep("TRUE", is.na(dat$V1))

      # If there are no clusters, return an empty data frame
      if(is.integer0(check) == TRUE) {
        sdat <- dat[, c(1, 4, 8, 11, 12)]
      } else {
        dat <- dat[-c(grep("TRUE", is.na(dat$V1))), ]
        sdat <- dat[, c(1, 4, 8, 11, 12)]
      }

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
