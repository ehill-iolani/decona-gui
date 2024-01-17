library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(plotly)
library(shinydashboard)
library(shinyalert)
library(DT)
library(ape)

# Accessory function to checks if an integer is of 0 length
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L
}

bc_sample <- function(x, y) {
  x <- x[sample(nrow(x), y, replace = FALSE), ]
  x <- x[!duplicated(x$gensp), ]
  x <- length(unique(x$gensp))
  return(x)
}

ui <- dashboardPage(skin = "red",
  dashboardHeader(title = "AIN Decona GUI"),
  dashboardSidebar(sidebarMenu(
    menuItem("Decona Classifier + Visualizer", tabName = "decona_class", icon = icon("upload")),
    menuItem("Decona Visualizer", tabName = "decona_viz", icon = icon("upload")),
    menuItem("Rarefaction Curve", tabName = "rarefaction", icon = icon("chart-line")),
    menuItem("Relative Abundance", tabName = "relab", icon = icon("chart-simple")),
    menuItem("BLAST Results", tabName = "blastres", icon = icon("table")),
    menuItem("Unclassified BLAST Hits", tabName = "unknown", icon = icon("question-circle")),
    menuItem("About", tabName = "about", icon = icon("info-circle"))
  )),
  # Formatting to make it iolani colors
  dashboardBody(tags$style(HTML("
  .box.box-solid.box-primary>.box-header {
    color:#fff;
    background:#000000
  }

  .box.box-solid.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .box.box-primary>.box-header {
    color:#000000;
    background:#fff
  }

  .box.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .skin-red .main-sidebar {
    background-color: #000000;
  }")),
  tabItems(
    tabItem(tabName = "decona_class",
            fluidRow(
            box(
              title = "Fastq Upload", status = "primary", solidHeader = TRUE, width = 12,
              fileInput("fastqs", "Input fastq files", accept = c("fastq", "fastq.gz"), multiple = TRUE),
              actionButton("fastqhelp", "Help!")),
            box(
              title = "Database Selection", status = "primary", solidHeader = TRUE, width = 12,
              selectInput("preloadedb", "Select a preloaded database (if you do not have one to upload)",
                          choices = c("I am providing my own database file", "Oahu Stream Fish", "Oahu Stream Decapod"),
                          selected = "I have provided my own"),
              fileInput("database", "Upload database .fasta file", accept = c("fasta")),
              actionButton("dbhelp", "Help!")),
            box(
              title = "Decona Parameters", status = "primary", solidHeader = TRUE, width = 12,
              numericInput("lowerlength", "Minimum amplicon length", value = 170),
              numericInput("upperlength", "Maximum amplicon length", value = 300),
              numericInput("quality", "Minimum Quality Score", value = 10),
              numericInput("clusterid", "Cluster Percent Identity", value = 0.95),
              numericInput("minclustersize", "Minimum Cluster Size", value = 10),
              numericInput("kmer", "Kmer-length", value = 10),
              numericInput("threads", "Threads", value = 2),
              actionButton("deconahelp", "Help!")),
            box(
              title = "Run Workflow", status = "primary", solidHeader = TRUE, width = 12,
              actionButton("run", "Run Decona Classifier and Visualizer!"),
              actionButton("dcvhelp", "Help!")))),
    tabItem(tabName = "decona_viz",
            fluidRow(
              box(
                title = "File Upload", status = "primary", solidHeader = TRUE, width = 12,
                fileInput("decona_results", "Upload Decona results file", accept = c("tsv")),
                actionButton("run_viz", "Run Decona Visualizer!"),
                actionButton("dvhelp", "Help!")))),
    tabItem(tabName = "rarefaction",
            fluidRow(
              box(
                title = "Rarefaction Curve",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("rarefaction"),
                downloadButton("downloadrarefaction", "Download Rarefaction Curve"),
                actionButton("rarefactionhelp", "Help!")))),
    tabItem(tabName = "relab",
            fluidRow(
              box(
                title = "Relative Abundance",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("relabres"),
                downloadButton("downloadfig", "Download Relative Abundance Plot"),
                actionButton("relabhelp", "Help!")))),
    tabItem(tabName = "blastres",
            fluidRow(
              box(
                title = "BLAST Results",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                DTOutput("blastres"),
                downloadButton("downloadblast", "Download BLAST Results"),
                actionButton("blastreshelp", "Help!")))),
    tabItem(tabName = "unknown",
            fluidRow(
              box(
                title = "Unclassified BLAST Hits",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                p("This tab will only be populated if the Decona Classifier + Visualizer was run. If the Decona Visualizer was run, re-upload the .fastq/fastq.gz files and run the Decona Classifier + Visualizer."),
                textOutput("unknownres"),
                downloadButton("downloadunknown", "Download Unclassified BLAST Hits"),
                actionButton("unknownhelp", "Help!")))),
    tabItem(tabName = "about",
            fluidRow(
              box(
                title = "About",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                p("This is a shiny app for the Decona pipeline.
                This app can used to classify amplicons sequenced by the Oxford Nanopore Technologies MinION platform or
                used to visualize clasifications previously generated by the Decona pipeline."),
                p("There are 2 primary workflows: Decona Classifier + Visualizer and Decona Visualizer."),
                p("The Decona Classifier + Visualizer workflow is used to classify and visualize amplicons sequenced by the Oxford Nanopore Technologies MinION platform.
                  The Decona Visualizer workflow is used to visualize classifications previously generated by the Decona pipeline.")))
  ))))

server <- shinyServer(function(input, output) {
  # Set large upload size limit (server side)q(0)
  options(shiny.maxRequestSize = 70 * 1024^2)

  ##########################################
  # EXECUTE DECONA CLASSIFIER + VISUALIZER #
  ##########################################
  observeEvent(input$run, {
    showModal(modalDialog("Running the Decona pipeline, please be patient...", footer = NULL))

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

    # Check if the files are fastq or fastq.gz; unzip if necessary
    setwd("/home/processing")
    if(any(str_detect(files, ".gz"))) {
      system("gunzip barcode*/*.gz")
      system("for i in */; do cd $i; for j in *; do mv $j ${j}.fastq; done; cd ..; done")
    }

    # Check the db selection and preps db accordingly
    if (input$preloadedb == "I am providing my own database file") {
      # Move the database file to the input directory
      file.copy(input$database$datapath, file.path("/home/processing"))
    } else if (input$preloadedb == "Oahu Stream Fish") {
      file.copy("/home/data/mifish_streamdb.fasta", "/home/processing/0.fasta")
    } else if (input$preloadedb == "Oahu Stream Decapod") {
      file.copy("/home/data/mideca_streamdb.fasta", "/home/processing/0.fasta")
    }

    # prep command and execute pipeline
    decona_command <- paste0("conda run -n decona decona -f -l ", input$lowerlength, " -m ", input$upperlength, " -q ", input$quality, " -c ", input$clusterid, " -n ", input$minclustersize, " -k ", input$kmer, " -T ", input$threads, " -B 0.fasta")
    system(decona_command, intern = TRUE)

    # Processes results
    # Detect how many barcodes were processed
    decona_out <- list.files("/home/processing/result/Racon")
    blast_out <- decona_out[grep("summary_BLAST", decona_out)]

    # Prep the data frame
    sdat <- data.frame()

    # Read in the raw text output
    for (i in blast_out) {
      dat <- read.delim(file.path("/home/processing/result/Racon", i), sep = "\t", header = FALSE, skip = 1, col.names = paste0("V",seq_len(100)), fill = TRUE)
      check <- grep("TRUE", is.na(dat$V1))

      # If there are no clusters, return an empty data frame
      if(is.integer0(check) == TRUE) {
        temp <- dat[, c(1, 4, 8, 11, 12)]
      } else {
        dat <- dat[-c(grep("TRUE", is.na(dat$V1))), ]
        temp <- dat[, c(1, 4, 8, 11, 12)]}

      # Rename the columns
      names(temp) <- c("reads", "percent_id", "e_value", "genus", "species")
      temp$percent_id[grep("TRUE", is.na(temp$percent_id))] <- 0

      # Add column for gensp and barcode
      temp$gensp <- paste(temp$genus, temp$species, sep = "_")
      temp$barcode <- str_extract(i, "barcode[0-9]+")

      # Add to the data frame
      sdat <- rbind(sdat, temp)
    }

    # Store sdat in output
    sdat <- as.data.frame(sdat)
    output$blastres <- renderDataTable(sdat)

    # Replace NA with unclassified
    sdat$gensp[is.na(sdat$e_value)] <- "unclassified"

    # Summarize reads per species
    condense <- sdat %>%
      group_by(barcode, gensp) %>%
      summarise(reads = sum(reads)) %>%
      arrange(desc(reads)) %>%
      mutate(gensp = factor(gensp, levels = unique(gensp)))

    # Calculate the relative abundance per unique barcode
    condense$rel_abund <- "fill"
    for (i in unique(condense$barcode)) {
      condense$rel_abund[condense$barcode == i] <- condense$reads[condense$barcode == i] / sum(condense$reads[condense$barcode == i])
    }

    # Make rel_abund numeric
    condense$rel_abund <- as.numeric(condense$rel_abund)

    # Display the relative abundance plot
    output$relabres <- renderPlotly({
      ggplotly(ggplot(condense, aes(x = barcode, y = rel_abund, fill = gensp)) +
                 geom_bar(stat = "identity") +
                 theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
                 labs(x = "Barcode", y = "Relative Abundance", fill = "Species"))
    })

    # Prepare figure for download
    output$downloadfig <- downloadHandler(
      filename = function() {
        paste("decona_relative_abundance_plot.pdf")
      },
      content = function(file) {
        pdf(file)
        print(ggplot(condense, aes(x = barcode, y = rel_abund, fill = gensp)) +
                 geom_bar(stat = "identity") +
                 theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
                 labs(x = "Barcode", y = "Relative Abundance", fill = "Species"))
      dev.off()
      }
    )
    # Prepare blast results for download
    output$downloadblast <- downloadHandler(
      filename = function() {
        paste("decona_blast_results.tsv")
      },
      content = function(file) {
        write.table(sdat, file, sep = "\t", row.names = FALSE)
      }
    )

    ###############################
    ### Find Unknown BLAST hits ###
    ###############################

    noid <- data.frame()

    # Read in the raw text output
    for (i in blast_out) {
      dat <- read.delim(file.path("/home/processing/result/Racon", i), sep = "\t", header = FALSE, skip = 1, col.names = paste0("V",seq_len(100)), fill = TRUE)
      check <- grep("TRUE", is.na(dat$V1))

      # If there are no clusters, return an empty data frame
      if(is.integer0(check) == TRUE) {
        temp <- dat[, c(1, 2, 4, 8, 11, 12)]
      } else {
        dat <- dat[-c(grep("TRUE", is.na(dat$V1))), ]
        temp <- dat[, c(1, 2, 4, 8, 11, 12)]
      }
      # Rename the columns
      names(temp) <- c("reads", "tag", "percent_id", "e_value", "genus", "species")
      temp$percent_id[grep("TRUE", is.na(temp$percent_id))] <- 0

      # Add columns for gensp, barcode
      temp$gensp <- paste(temp$genus, temp$species, sep = "_")
      temp$barcode <- str_extract(i, "barcode[0-9]+")

      # Add to the data frame
      noid <- rbind(noid, temp)
    }

    # Extracts unclassified reads from noid
    noid <- as.data.frame(noid)
    noid$gensp[is.na(noid$e_value)] <- "unclassified"
    noid <- noid[noid$gensp == "unclassified", ]

    # Add catch statment if there are no unclassified BLAST hits
    if (nrow(noid) == 0) {
      output$unknownres <- renderText({
        paste("No clusters returned an unclassified status!")
      })
      removeModal()
    } else {
      # Reconstuct file name
      noid$filename <- paste0("polished-", noid$reads, "-", noid$tag, ".fasta")

      # Set the working directory to the data directory
      # Cannot create empty DNAbin object, so need to read in some
      # sample data to get it to work
      setwd("/home/processing/data")
      dirs <- list.dirs()[str_detect(list.dirs(), "barcode[0-9]+_concatenated/multi-seq")]
      data(woodmouse)
     fastas <- as.list(woodmouse)
      for (i in dirs) {
        fasta_files <- list.files(i, pattern = "n-polished.*.fasta")
        for (j in fasta_files) {
          temp <- read.FASTA(file.path(i, j))
          fastas <- c(fastas, temp)
        }
      }
      fastas <- fastas[-c(1:15)]
      # names(fastas) <- str_replace_all(names(fastas), "polished-", "polished_")
      fastas <- fastas[names(fastas) %in% noid$filename]

      # add the matching barcode from noid to the names(fastas) based on noid$filename
      for (i in 1:length(names(fastas))) { # nolint
        index <- noid$filename %in% names(fastas)[i]
        name <- paste(noid$barcode[index], names(fastas)[i], sep = "_")
        names(fastas)[i] <- name
      }

      # Show the unknown BLAST hits
      output$unknownres <- renderText({
        paste("The following unknown BLAST hits were found, please download and manually BLAST:", paste(names(fastas), collapse = ", "))
      })

      # Download for unknown BLAST hits
      output$downloadunknown <- downloadHandler(
        filename = function() {
          paste("decona_unknown_blast_hits.fasta")
        },
        content = function(file) {
          write.FASTA(fastas, file)
        })
    }

    # Rarefaction analysis + plot
    # make a vector containing all the reads per unique species by barcode
    condenser <- as.data.frame(condense)
    condenser$gensp <- as.character(condenser$gensp)
    rare <- c()
    for (i in 1:length(condense$barcode)) { # nolint
      temp <- rep(condense$gensp[i], condense$reads[i])
      temp <- data.frame(gensp = temp)
      temp$barcode <- condense$barcode[i]
      rare <- rbind(rare, temp)
    }

    # Randomly sample 1:length(reads) and calculate the number of unique species
    # for each sample
    rarefaction <- data.frame()
    for (i in unique(rare$barcode)) {
      temp <- rare[rare$barcode == i, ]
      for (j in seq(1, nrow(temp), by = 10)) {
        temp2 <- mean(replicate(2, bc_sample(temp, j)))
        temp2 <- as.data.frame(temp2)
        temp2$reads <- j
        temp2$barcode <- i
        rarefaction <- rbind(rarefaction, temp2)
      }
    }
    names(rarefaction) <- c("unique_species", "reads", "barcode")

    # Plot the rarefaction curve by barcode
    output$rarefaction <- renderPlotly({
      ggplotly(ggplot(rarefaction, aes(x = reads, y = unique_species, color = barcode)) +
                 geom_point() +
                 geom_line() +
                 labs(x = "Number of Reads", y = "Number of Unique Species", color = "Barcode"))
    })

    # Prepare rarefaction curve for download
    output$downloadrarefaction <- downloadHandler(
      filename = function() {
        paste("decona_rarefaction_curve.pdf")
      },
      content = function(file) {
        pdf(file)
        print(ggplot(rarefaction, aes(x = reads, y = unique_species, color = barcode)) +
                 geom_point() +
                 geom_line() +
                 labs(x = "Number of Reads", y = "Number of Unique Species", color = "Barcode"))
        dev.off()
      }
    )

    # Indicate that the process is complete
    removeModal()
    shinyalert(
      title = "Completed!",
      text = "Decona Classifier + Visualizer has finished running! Results can be viewed under the Relative Abundance and BLAST results tabs.",
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
    )
  })

  #############################
  # EXECUTE DECONA VISUALIZER #
  #############################
  observeEvent(input$run_viz, {
    showModal(modalDialog("Running the Decona Visualizer, please be patient...", footer = NULL))

    # Read in the results file
    decona_results <- read.delim(input$decona_results$datapath, sep = "\t")
    sdat <- decona_results

    # Display the blast results
    output$blastres <- renderDataTable(sdat)

    # Summarize reads per species
    condense <- sdat %>%
      group_by(barcode, gensp) %>%
      summarise(reads = sum(reads)) %>%
      arrange(desc(reads)) %>%
      mutate(gensp = factor(gensp, levels = unique(gensp)))

    # Calculate the relative abundance per unique barcode
    condense$rel_abund <- "fill"
    for (i in unique(condense$barcode)) {
      condense$rel_abund[condense$barcode == i] <- condense$reads[condense$barcode == i] / sum(condense$reads[condense$barcode == i])
    }

    # Make rel_abund numeric
    condense$rel_abund <- as.numeric(condense$rel_abund)

    # Display the relative abundance plot
    output$relabres <- renderPlotly({
      ggplotly(ggplot(condense, aes(x = barcode, y = rel_abund, fill = gensp)) +
                 geom_bar(stat = "identity") +
                 theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
                 labs(x = "Barcode", y = "Relative Abundance", fill = "Species"))
    })

    # Prepare figure for download
    output$downloadfig <- downloadHandler(
      filename = function() {
        paste("decona_relative_abundance_plot.pdf")
      },
      content = function(file) {
        pdf(file)
        print(ggplot(condense, aes(x = barcode, y = rel_abund, fill = gensp)) +
                 geom_bar(stat = "identity") +
                 theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
                 labs(x = "Barcode", y = "Relative Abundance", fill = "Species"))
        dev.off()
      }
    )
    # Prepare blast results for download
    output$downloadblast <- downloadHandler(
      filename = function() {
        paste("decona_blast_results.tsv")
      },
      content = function(file) {
        write.table(sdat, file, sep = "\t", row.names = FALSE)
      }
    )

    # Indicate that the process is complete
    removeModal()
    shinyalert(
      title = "Completed!",
      text = "Decona Visualizer has finished running! Results can be viewed under the Relative Abundance and BLAST results tabs.",
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
    )
  })

  #################
  # HELP MESSAGES #
  #################
  observeEvent(input$fastqhelp, {
    shinyalert(
      title = "Fastq Upload Help",
      html = TRUE,
      text = "<b>Please upload your .fastq/.fastq.gz files from your ONT MinION run and a database .fasta file.</b><br><br>
      Following the sequencing run, make sure to collected all of the relevant .fastq/.fastq.gz files from the barcode directories.<br><br>
      Once you have collected all of the relevant .fastq/.fastq.gz files, upload them here.<br><br>",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })

  observeEvent(input$dbhelp, {
    shinyalert(
      title = "Database Selection Help",
      html = TRUE,
      text = "<b>The database .fasta file should contain the gene sequences of the species you believe are present in your sample site(s).</b><br><br>
      If you do not have a database .fasta file, you can select one of the preloaded databases.<br><br>
      Make sure to select the correct database for your amplicon type.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })

  observeEvent(input$deconahelp, {
    shinyalert(
      title = "Decona Parameters Help",
      html = TRUE,
      text = "For most workflows, the default parameters should work, however, <b>make sure to adjust the upper and lower limits of the amplicon length accordingly.</b><br><br>
      Set the minimum and maximum amplicon lengths, minimum quality score, cluster percent identity, minimum cluster size, kmer-length, and number of threads.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })

  observeEvent(input$dcvhelp, {
    shinyalert(
      title = "Decona Classifier + Visualizer Help",
      html = TRUE,
      text = "To run the Decona Classifer + Visualizer, <b>please upload your .fastq/.fastq.gz files from your ONT MinION run and a database .fasta file.</b><br><br>
      The database .fasta file should contain the gene sequences of the species you believe are present in your sample site(s).<br><br>
      Set the minimum and maximum amplicon lengths, minimum quality score, cluster percent identity, minimum cluster size, kmer-length, and number of threads.<br><br>
      <b>The Decona Classifier + Visualizer will then classify your reads and display the results under the Rarefaction Curve, Relative Abundance and BLAST results tabs.</b>",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
  observeEvent(input$dvhelp, {
    shinyalert(
      title = "Decona Visualizer Help",
      html = TRUE,
      text = "To run the Decona Visualizer, please <b>upload your Decona results file from a previous Decona Classifer + Visualizer run.</b><br><br>
      <b>The Decona Visualizer will then display the results under the Rarefaction Curve, Relative Abundance and BLAST results tabs.</b>",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
  observeEvent(input$rarefactionhelp, {
    shinyalert(
      title = "Rarefaction Curve Help",
      html = TRUE,
      text = "<b>The Rarefaction Curve tab displays the number of unique species found in your sample(s) as a function of the number of reads.</b><br><br>
      The Rarefaction Curve tab is populated by the Decona Classifier + Visualizer and Decona Visualizer workflows.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
  observeEvent(input$relabhelp, {
    shinyalert(
      title = "Relative Abundance Help",
      html = TRUE,
      text = "<b>The Relative Abundance tab displays the relative abundance of each species found in your sample(s).</b><br><br>
      The Relative Abundance tab is populated by the Decona Classifier + Visualizer and Decona Visualizer workflows.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
  observeEvent(input$blastreshelp, {
    shinyalert(
      title = "BLAST Results Help",
      html = TRUE,
      text = "<b>The BLAST Results tab displays the BLAST results for each species found in your sample(s).</b><br><br>
      The BLAST Results tab is populated by the Decona Classifier + Visualizer and Decona Visualizer workflows.<br><br>
      <b>You can download these results by clicking the Download BLAST Results button.</b><br><br>
      You can then upload the downloaded results and reanalyze them in the future using the Decona Visualizer workflow.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
  observeEvent(input$unknownhelp, {
    shinyalert(
      title = "Unknown BLAST Hits Help",
      html = TRUE,
      text = "<b>The Unclassified BLAST Hits tab displays the unclassified BLAST hits found in your sample(s).</b><br><br>
      The Unclassified BLAST Hits tab is populated by the Decona Classifier + Visualizer workflow.<br><br>
      <b>You can download these unclassified sequences and manually BLAST them to determine their identity by clicking the Download Unclassified BLAST Hits button.</b>",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#000000"
    )
  })
})

shinyApp(ui, server)
