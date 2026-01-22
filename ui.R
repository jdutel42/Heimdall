###################################################################################"
# ---- Initialization ---- #

#############
# Libraries #
#############

library(shiny)
library(ggplot2)
library(SingleCellExperiment)
library(shinycssloaders)

###########
# Options #
###########

options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

###################################################################################"
# ---- User Interface ---- #

ui <- fluidPage(
  
  
  ###################
  # ---- Title ---- #
  ###################
  
  titlePanel("Single-cellExperiment FeaturePlot Explorer"),
  
  
  #####################
  # ---- Sidebar ---- #
  #####################
  
  sidebarLayout(
    sidebarPanel(
      
      ########################
      # ---- Input file ---- #
      ########################
      
      fileInput(
        "sce_rds",
        "Upload SCE object (.rds)",
        accept = ".rds"
      ),
      
      ###################
      # ---- Assay ---- #
      ###################
      
      shinycssloaders::withSpinner(
        uiOutput("assay_ui")
      ),
      
      ###################
      # ---- Genes ---- #
      ###################
      
      shinycssloaders::withSpinner(
        uiOutput("gene_ui")
      ),
      
      #######################
      # ---- Reduction ---- #
      #######################
      
      selectInput(
        "reduction",
        "Embedding",
        # choices = c("UMAP_uwot", "TSNE", "PCA"),
        choices = c("PCA", "TSNE", "UMAP"),
        selected = "UMAP"
      ),
      
      
      #########################
      # ---- Plot Tuning ---- #
      #########################
      
      sliderInput(
        "pt_size",
        "Point size",
        min = 0.1,
        max = 2,
        value = 0.6,
        step = 0.1
      ),
      
      
      ##################################
      # ---- Plot Tuning Metadata ---- #
      ##################################
      
      shinycssloaders::withSpinner(
        uiOutput("metadata_ui")
        ),
    
      
      
      
      # selectInput(
      #   "metadata",
      #   "Choose a metadata to display",
      #   nested_list_metadata
      # ),
      
      
      ###########################
      # ---- Download Plot ---- #
      ###########################
      
      # Button to download the PCA plot
      downloadButton(
        "download_plot", 
        "Download Plot"
      )
      
      
      
      
      
    ),
    
    ################################
    # ---- Plot visualization ---- #
    ################################
    
    mainPanel(
      shinycssloaders::withSpinner(
        plotOutput(
          "featureplot",
          height = "750px"
        )
      )
    )
  )
)
