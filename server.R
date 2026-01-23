server <- function(input, output, session) {

  
  #############################
  # ---- Load SCE object ---- #
  #############################
  
  sce_obj <- reactive({
    
    req(input$sce_rds)
    
    sce <- readRDS(input$sce_rds$datapath)
    
    validate(
      need(
        inherits(sce, "SingleCellExperiment"),
        "The uploaded file is not a SingleCellExperiment object"
      )
    )
    
    sce
  })


  #############################
  # ---- Assay selection ---- #
  #############################
  
  output$assay_ui <- renderUI({
    
    req(input$sce_rds)
    
    selectInput(
      "assay",
      "2. Which assay ?",
      choices = assayNames(sce_obj()),
      selected = assayNames(sce_obj())[1]
    )
  })
  
  
  #################################
  # ---- Embedding selection ---- #
  #################################
  
  output$embedding_ui <- renderUI({
    
    req(input$assay)
    
    selectInput(
      "embedding",
      "3. Which embedding ?",
      choices = reducedDimNames(sce_obj()),
      selected = reducedDimNames(sce_obj())[1]
    )
  })

  ###############################
  # ---- Feature selection ---- #
  ###############################
  
  output$feature_ui <- renderUI({
    
    req(input$embedding)
    
    selectInput(
      "feature",
      "4. Which features ?",
      # choices = c("UMAP_uwot", "TSNE", "PCA"),
      choices = c("Expression", "Clusters", "Others Metadata"),
      selected = "Expression"
    )
  })


  ############################
  # ---- Gene selection ---- #
  ############################

  output$gene_ui <- renderUI({
    req(sce_obj(), input$embedding)
    genes <- rownames(assay(sce_obj(), input$assay))
    selectizeInput(
      "genes",
      "Which Gene(s) ?",
      choices = genes,
      multiple = TRUE,
      options = list(
        maxItems = 4,
        placeholder = "Type gene name(s)..."
      )
    )
  })
  
  
  ######################
  # ---- Metadata ---- #
  ######################
  
  # # Create a nested list to store all possible metavalues in each col
  output$metadata_ui <- renderUI({
    req(sce_obj())
    
    # Store names of metadata present in the SCE object
    colnames_metadata_list <- colnames(colData(sce_obj()))
    # If col have more than 20 unique value are present, remove them from the choices
    # colnames_metadata_list <- colnames_metadata_list[sapply(colnames_metadata_list, function(col) {
    #   length(unique(colData(sce_obj())[[col]])) <= 30
    # })]
    # colnames_metadata_list <- setdiff(colnames_metadata_list, c("nCounts", "nFeatures"))
    
    selectInput(
      "metadata_col",
      "Choose data to color by",
      choices = colnames_metadata_list
    )
  })
  
  # output$metadata_val_ui <- renderUI({
  #   req(sce_obj(), input$metadata_col)
  #   
  #   vals <- unique(colData(sce_obj())[[input$metadata_col]])
  #   
  #   selectInput(
  #     "metadata_val",
  #     "Choose metadata value",
  #     choices = vals
  #   )
  # })
  
    # selectInput(
    #   "metadata_val",
    #   "Choose a metadata to display",
    #   choices = setNames(
    #     lapply(colnames_metadata_list, function(col) {
    #       unique(colData(sce_obj())[[col]])
    #     }),
    #     colnames_metadata_list
    #   ),
    # )



  
  #########################
  # ---- FeaturePlot ---- #
  #########################
  
  ## Calculate
  featureplot_obj <- reactive({
    sce <- sce_obj()
    validate(
      need(
        input$embedding %in% reducedDimNames(sce),
        paste(
          "embedding",
          input$embedding,
          "not found in object. Ask bioinfo to compute it."
        )
      )
    )
    
    # Create the data.frame for ggplot
    
    ## Coordinates embedding
    emb <- as.data.frame(reducedDim(sce, input$embedding))
    colnames(emb)[1:2] <- c("Dim1", "Dim2")
    
    ## Expression values
    expr_assay <- assay(sce, input$assay)[input$genes, , drop = FALSE]
    ### Convert sparse matrix -> dense matrix
    expr_mat <- as.matrix(expr_assay)
    expr <- as.data.frame(t(expr_mat))
    
    ## Metadata
    meta <- as.data.frame(colData(sce))
    
    ## Combine
    df <- cbind(
      emb,
      expr,
      meta
    )
    
    df_long <- tidyr::pivot_longer(
      df,
      cols = all_of(input$genes),
      names_to = "gene",
      values_to = "Expression"
    )
    
    #   ggplot2::ggplot(
    #     df_long,
    #     ggplot2::aes(x = Dim1, y = Dim2, color = Expression)
    #   ) +
    #     ggplot2::geom_point(size = input$pt_size) +
    #     ggplot2::scale_color_gradient(low = "lightgrey", high = "blue") +
    #     ggplot2::facet_wrap(~ gene) +
    #     ggplot2::theme_classic(base_size = 16, base_line_size = 0.5, base_rect_size = 0.5, ink = "black")
    #     # ggplot2::theme_bw()
    # })
    
    
    
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = Dim1, y = Dim2)
    ) +
      ggplot2::geom_point(size = input$pt_size) +
      ggplot2::facet_wrap(~ gene) +
      ggplot2::theme_classic(
        base_size = 16,
        base_line_size = 0.5,
        base_rect_size = 0.5,
        ink = "black"
      )
    
    print(df_long)
    
    if (is.null(input$metadata_col) || input$metadata_col == "") {
      
      ## Par défaut : expression
      p <- p +
        ggplot2::aes(color = Expression) +
        ggplot2::scale_color_gradient(
          low = "lightgrey",
          high = "blue"
        )
      
    } else {
      
      ## Metadata sélectionnée
      p <- p +
        ggplot2::aes(color = .data[[input$metadata_col]])
      
      # ## Optionnel : adapter l’échelle selon le type
      # if (is.numeric(df_long[[input$metadata_col]])) {
      #   p <- p + ggplot2::scale_color_viridis_c()
      # } else {
      #   p <- p + ggplot2::scale_color_brewer(palette = "Set2")
      # }
    }
    
    p
  })
  

  ## Ploting featureplot
  output$featureplot <- renderPlot({
    req(sce_obj(), input$genes, input$assay, input$embedding)
    featureplot_obj()
  })
  

  ###########################
  # ---- Download plot ---- #
  ###########################
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("FeaturePlot_SCE_", Sys.Date(), "_", format(Sys.time(), "%X"), ".png", sep = "") # To improve to dynamically change the name (eg gene name or whatever)
    },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 150)
      print(featureplot_obj())
      dev.off()
    }
  )


  
  
  
  
  
  
  
  
  # nested_list_metadata <- reactive({
  #   
  #   # Store names of metadata present in the SCE object
  #   colnames_metadata_list <- colnames(colData(sce_obj()))
  # 
  #   # Create the list of list with names according colnames
  #   setNames(
  #     lapply(colnames_metadata_list, function(col) {
  #       unique(colData(sce_obj())[[col]])
  #     }),
  #     colnames_metadata_list
  #   )
  # })
  
  
  
}
