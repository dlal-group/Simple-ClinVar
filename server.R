library(shiny) 
library(ggplot2) 
library(readr) 
library(stringr)
library(tidyr)
library(Hmisc)
library(ggrepel)
library(DT)
library(ggrepel)
library(dplyr)
library(tm)
library(shinycssloaders)
library("shinythemes")
# With Drawprotein
library(drawProteins)
library(magrittr)
#library(rlang)

# Load main data ----------------------------------------------------------
clean     <- read_csv("clinvar.feb19_2.pf")
genes     <- read_delim("genes.refSeq", "\t", escape_double = FALSE, trim_ws = TRUE)
sequences <- read_delim("gene-ccds-seq-length-uniprot.txt", "\t",
                        escape_double = FALSE,
                        col_types = cols(Length = col_number()),
                        trim_ws = TRUE)

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

shinyServer(function(input, output, session) {
  
  # Cambios de pestagna
  # Principal 
  observeEvent(input$goButton,  { updateTabsetPanel(session, "inTabset", selected = "panel2")})
  observeEvent(input$newsearch, { updateTabsetPanel(session, "inTabset", selected = "panel1")})
  # Resultados
  observeEvent(input$seevariants,    { updateTabsetPanel(session, "displaywhat", selected = "1")})
  observeEvent(input$seegenes,       { updateTabsetPanel(session, "displaywhat", selected = "2")})
  observeEvent(input$seephenotypes,  { updateTabsetPanel(session, "displaywhat", selected = "3")})
  observeEvent(input$seetable,       { updateTabsetPanel(session, "displaywhat", selected = "4")})  
  
  query <- eventReactive(input$goButton,{
      return(toupper(input$genename))  
  })
   
    # poner busqueda en mayuscula
  output$displayquery <- renderText({
    paste0("\"", query(), "\"")
  })
  
  subset_one <- reactive({
      if (query() == "CLINVAR") {
      clean
    } else if (query() %in% genes$gene) {
      subset(clean, grepl(query(), clean$GeneSymbol, ignore.case = TRUE))
    } else {
      validate(need(grepl(query(), clean$PhenotypeList, ignore.case=TRUE), "Please repeat your search."))
      subset(clean, grepl(query(), clean$PhenotypeList, ignore.case=TRUE))
    }
  }) # es gen o no. primer grepeo
  
  # Dinamic filtering generates subset_two variable -------------------------------------------------------
  
  clinicalinput    <- reactive({input$clinicalinput})
  consequenceinput <- reactive({input$consequenceinput})
  typeinput <- reactive({input$typeinput})
  reviewinput <- reactive({input$reviewinput})
  
  query_two <- eventReactive(input$filter|input$goButton,{
    return(c(clinicalinput(),consequenceinput(),typeinput(),reviewinput()))
  })
  
  clinical_index <- eventReactive(input$goButton|input$filter,{
    clin_index <- c()
    if (length(query_two()) == 0) {
      clin_index <- seq(1,nrow(subset_one()))}
    else if (length(clinicalinput())== 0){
      clin_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in clinicalinput()){
        if (i == "Pathogenic"){
          clin_index <- c(clin_index, grep("6",subset_one()$ClinicalSignificance_grouped))}
        if (i == "Likely pathogenic"){
          clin_index <- c(clin_index, grep("5",subset_one()$ClinicalSignificance_grouped))}
        if (i == "Risk factor/Association"){
          clin_index <- c(clin_index, grep("4",subset_one()$ClinicalSignificance_grouped))}
        if (i == "Uncertain/conflicting"){
          clin_index <- c(clin_index, grep("3",subset_one()$ClinicalSignificance_grouped))}
        if (i == "Protective/Likely benign"){
          clin_index <- c(clin_index, grep("2",subset_one()$ClinicalSignificance_grouped))}
        if (i == "Benign"){
          clin_index <- c(clin_index, grep("1",subset_one()$ClinicalSignificance_grouped))}
        }}
    return(clin_index)
  })
  
  consequence_index <- eventReactive(input$goButton|input$filter,{
    con_index <- c()
    if (length(query_two()) == 0){
      con_index <- seq(1,nrow(subset_one()))}
    else if (length(consequenceinput())== 0){
      con_index <- seq(1,nrow(subset_one()))}
    else {
    for (i in consequenceinput()){
      if (i == "Missense"){
        con_index <- c(con_index, grep("Missense", subset_one()$consequence))}
      if (i == "Stop-gain"){
        con_index <- c(con_index, grep("Stop gain", subset_one()$consequence))}
      if (i == "In frame indel"){
        con_index <- c(con_index, grep("In frame indel", subset_one()$consequence))}
      if (i == "Frameshift"){
        con_index <- c(con_index, grep("Frameshift", subset_one()$consequence))}
      if (i == "Synonymous"){
        con_index <- c(con_index, grep("Synonymous", subset_one()$consequence))}
        }}
    return(con_index)
  })
  
  type_index <- eventReactive(input$goButton|input$filter,{
    type_index <- c()
    if (length(query_two()) == 0){
      type_index <- seq(1,nrow(subset_one()))}
    else if (length(typeinput())== 0){
      type_index <- seq(1,nrow(subset_one()))}
    else {
    for (i in typeinput()){
      if (i == "SNV"){
        type_index <- c(type_index, grep("SNV", subset_one()$Type))}
      if (i == "Deletion"){
        type_index <- c(type_index, grep("Deletion", subset_one()$Type))}
      if (i == "Duplication"){
        type_index <- c(type_index, grep("Duplication", subset_one()$Type))}
      if (i == "Deletion"){
        type_index <- c(type_index, grep("Deletion", subset_one()$Type))}
      if (i == "Indel"){
        type_index <- c(type_index, grep("Indel", subset_one()$Type))}
      if (i == "Inversion"){
        type_index <- c(type_index, grep("Inversion", subset_one()$Type))}
        }}
    return(type_index)
  })
  
  
  review_index <- eventReactive(input$goButton|input$filter,{
    review_index <- c()
    if (length(query_two()) == 0){
      review_index <- seq(1,nrow(subset_one()))}
    else if (length(reviewinput())== 0){
      review_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in reviewinput()){
        if (i == "No assertion criteria provided"){
          review_index <- c(review_index, grep("No assertion criteria provided", subset_one()$review))}
        if (i == "Criteria provided/ conflicting interpretations"){
          review_index <- c(review_index, grep("Criteria provided/ conflicting interpretations", subset_one()$review))}
        if (i == "Criteria provided/ single submitter"){
          review_index <- c(review_index, grep("Criteria provided/ single submitter", subset_one()$review))}
        if (i == "Practice guideline"){
          review_index <- c(review_index, grep("Practice guideline", subset_one()$review))}
        if (i == "Reviewed by expert panel"){
          review_index <- c(review_index, grep("Reviewed by expert panel", subset_one()$review))}
        if (i == "Criteria provided/ multiple submitters/ no conflicts"){
          review_index <- c(review_index, grep("Criteria provided/ multiple submitters/ no conflicts", subset_one()$review))}
      }}
    return(review_index)
  })
  
  merged_index <- reactive({
    sub_index <- c()
      sub_index <- intersect(intersect(intersect(consequence_index(),clinical_index()),type_index()),review_index())
    return(sub_index)
  })
  
  subset_two <- eventReactive(input$goButton|input$filter,{
    if (length(merged_index()) == 0) {
      h <- NULL
      x <- as.data.frame(h)
      return(x)
    } else {
      x <- subset_one()[merged_index(),]
      return(x)}
    })
  
  # N of observations -------------------------------------------------------
  numberobservations <- reactive({nrow(subset_one())})
  
  output$displayobservation <- renderText({
    if (nrow(subset_two()) == 0) {
      return(paste0("0 variants"))
    } else {
      return(paste0(nrow(subset_two())," variants"))
    }
  })
  
  
  # histogram "review" ----------------------------------------------------------
  allreview <-   reactive({unique(clean$review)})
  allreview_1 <- reactive({paste0(allreview(), collapse=',')})
  allreview_2 <- reactive({str_split_fixed(allreview_1(), ',', Inf)})
  allreview_3 <- reactive({as.data.frame(t(allreview_2()))})
  allreview_4 <- reactive({unique(allreview_3())})
  allreview_5 <- reactive({as.character(allreview_4()$V1)})
  #### JOINING ALL REVIEW
  reviewgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$review, collapse=',')})
  reviewgene_2 <- reactive({str_split_fixed(reviewgene(), ',', Inf)})
  reviewgene_3 <- reactive({as.data.frame(t(reviewgene_2()))})
  #### SORTING LEVELS AND ADDING LEVELS SO THAT ALL LEVELS ARE ALWAYS IN EACH PLOT
  reviewgene_4 <- reactive({factor(reviewgene_3()$V1, levels = allreview_5())})
  reviewgene_5 <- reactive({data.frame(reviewgene_4())})
  reviewgene_6 <- reactive({table(reviewgene_5())})
  reviewgene_7 <- reactive({sort(reviewgene_6())})
  reviewgene_8 <- reactive({names(reviewgene_7())})
  reviewgene_9 <- reactive({factor(reviewgene_5()$reviewgene_4.., levels = reviewgene_8())})
  reviewgene_10 <- reactive({data.frame(reviewgene_9())}) 
  ####GGPLOT REVIEW
  output$review <- renderPlot({
      if (nrow(subset_two())==0){
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.34, y = 0.9, paste(""), 
             cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
      ggplot(reviewgene_10(),mapping = aes(x=reviewgene_10()$reviewgene_9..)) +
        ylab("Counts") + 
        xlab("") +
        geom_bar(fill = "#43AC6A", show.legend = F) +
        scale_x_discrete(drop=FALSE, labels=c("No assertion criteria provided"="No assertion\ncriteria provided", 
                                              "no assertion criteria provided"="no assertion\ncriteria provided",
                                              "Criteria provided/ conflicting interpretations"= "Criteria provided/\nconflicting\ninterpretations", 
                                              "Criteria provided/ single submitter"="Criteria provided/\nsingle submitter", 
                                              "Practice guideline"="Practice\nguideline",
                                              "Reviewed by expert panel"="Reviewed by\nexpert panel",
                                              "Criteria provided/ multiple submitters/ no conflicts"="Criteria provided/\nmultiple submitters/\nno conflicts"))+
        scale_y_continuous(labels = scales::comma)+
        coord_flip()+
        theme_light()+
        theme(title = element_blank())
    }
  })
  
  # histogram type ----------------------------------------------------------
  alltype <-   reactive({unique(clean$Type)})
  alltype_1 <- reactive({paste0(alltype(), collapse=',')})
  alltype_2 <- reactive({str_split_fixed(alltype_1(), ',', Inf)})
  alltype_3 <- reactive({as.data.frame(t(alltype_2()))})
  alltype_4 <- reactive({unique(alltype_3())})
  alltype_5 <- reactive({as.character(alltype_4()$V1)})
  #### JOINING ALL TYPE ANNOTATIONS (FROM OF OBSERVATIONS FOR ONE GENE/DISEASE)
  
  typegene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$Type, collapse=',')})
  typegene_2 <- reactive({str_split_fixed(typegene(), ',', Inf)})
  typegene_3 <- reactive({as.data.frame(t(typegene_2()))})
  #### SORTING LEVELS AND ADDING LEVELS SO THAT ALL LEVELS ARE ALWAYS IN EACH PLOT
  typegene_4 <- reactive({factor(typegene_3()$V1, levels = alltype_5())})
  typegene_5 <- reactive({data.frame(typegene_4())})
  typegene_6 <- reactive({table(typegene_5())})
  typegene_7 <- reactive({sort(typegene_6())})
  typegene_8 <- reactive({names(typegene_7())})
  typegene_9 <- reactive({factor(typegene_5()$typegene_4.., levels = typegene_8())})
  typegene_10 <- reactive({data.frame(typegene_9())}) 
  ####GGPLOT TYPE
  output$Type <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste("No variants 
                    remaining after filtering."), 
           cex = 1.5, col = "black", family="sans", font=1, adj=0.5)
    } else {
      ggplot(typegene_10(),mapping = aes(x=typegene_10()$typegene_9..)) +
        ylab("Counts") + 
        xlab("") +
        geom_bar(fill = "#43AC6A", show.legend = F) +
        scale_x_discrete(drop=FALSE, labels=c("Undetermined variant"="Undetermined\nvariant"))+
        scale_y_continuous(labels = scales::comma)+
        coord_flip()+
        theme_light()+
        theme(title = element_blank())
    }
  })
  
  # histogram consequence ---------------------------------------------------
  #### ALL POSSIBLE CONSEQUENCE ANNOTATIONS FROM THE BIG DATASET (CLEAN)-> LEVELS FOR PLOT: allcons_5
  allcons   <- reactive({unique(clean$consequence)})
  allcons_1 <- reactive({paste0(allcons(), collapse=',')})
  allcons_2 <- reactive({str_split_fixed(allcons_1(), ',', Inf)})
  allcons_3 <- reactive({as.data.frame(t(allcons_2()))})
  allcons_4 <- reactive({unique(allcons_3())})
  allcons_5 <- reactive({as.character(allcons_4()$V1)})
  
  #### JOINING ALL CONSEQUENCE ANNOTATIONS  (FROM ALL OBSERVATIONS OF ONE GENE/DISEASE)
  consgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$consequence, collapse=',')})    
  consgene_2 <- reactive({str_split_fixed(consgene(), ',', Inf)})
  consgene_3 <- reactive({as.data.frame(t(consgene_2()))})
  
  #### SORTING LEVELS AND ADDING LEVELS SO THAT ALL LEVELS ARE ALWAYS IN EACH PLOT
  consgene_4 <- reactive({factor(consgene_3()$V1, levels = allcons_5())})
  consgene_5 <- reactive({data.frame(consgene_4())})
  consgene_6 <- reactive({table(consgene_5())})
  consgene_7 <- reactive({sort(consgene_6())})
  consgene_8 <- reactive({names(consgene_7())})
  consgene_9 <- reactive({factor(consgene_5()$consgene_4.., levels = consgene_8())})
  consgene_10 <- reactive({data.frame(consgene_9())}) 
  
  ### GGPLOT CONSEQUENCE
  output$consequence <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste(""), 
           cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
    ggplot(consgene_10(),mapping = aes(x=consgene_10()$consgene_9..)) +
      ylab("Counts") + 
      xlab("") +
      geom_bar(fill = "#43AC6A", show.legend = F)+
      scale_x_discrete(drop=FALSE)+
      scale_y_continuous(labels = scales::comma)+
      coord_flip() +
      theme_light() +
      theme(title = element_blank())
    }
  })
  # histogram clinical significance -----------------------------------------
  #### ALL POSSIBLE CLINICAL SIGNFICANCE ANNOTATIONS FROM THE BIG DATASET-> LEVELS FOR PLOT: allclinicals_6
  allclinicals_6 <- reactive({c("1","2","3","4","5","6","NA")}) 
  
  #### JOINING ALL CLINICAL SIGNIFICANCE ANNOTATIONS FROM ONE GENE: clinicalgene_3$V1
  clinicalgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$ClinicalSignificance_grouped, collapse=',')})
  clinicalgene_2 <- reactive({str_split_fixed(clinicalgene(), ',', Inf)})
  clinicalgene_3 <- reactive({as.data.frame(t(clinicalgene_2()))})

  
  #### SORTING LEVELS AND ADDING LEVELS SO THAT ALL LEVELS ARE ALWAYS IN EACH PLOT
  clinicalgene_4 <- reactive({factor(clinicalgene_3()$V1, levels = allclinicals_6())})
  clinicalgene_5 <- reactive({data.frame(clinicalgene_4())})
  clinicalgene_6 <- reactive({table(clinicalgene_5())})
  clinicalgene_7 <- reactive({sort(clinicalgene_6())})
  clinicalgene_8 <- reactive({names(clinicalgene_7())})
  clinicalgene_9 <- reactive({factor(clinicalgene_5()$clinicalgene_4.., levels = clinicalgene_8())})
  clinicalgene_10 <- reactive({data.frame(clinicalgene_9())}) 
  clinicalgene_10a <- reactive({clinicalgene_10()$clinicalgene_9[is.na(clinicalgene_10()$clinicalgene_9) == F]}) 
  clinicalgene_10b  <- reactive({data.frame(clinicalgene_10a())}) 
  #### GGPLOT CLINICAL SIGNIFICANCE
  output$ClinicalSignificance <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste(""), 
           cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
    ggplot(clinicalgene_10b(), aes(x=clinicalgene_10b()$clinicalgene_10a)) +
      ylab("Counts") + 
      xlab("") +
      geom_bar(fill = "#43AC6A", show.legend = F) +
      scale_x_discrete(drop=FALSE, labels=c("1"="Benign", "2"="Protective/\nLikely benign", "3"="Uncertain/\nconflicting", 
                                            "4"="Risk factor/\nAssociation","5"="Likely\npathogenic","6"="Pathogenic","NA"="NA"))+
      scale_y_continuous(labels = scales::comma)+
      coord_flip() +
      theme_light() +
      theme(title = element_blank())
    }
  })
  
  # Protein Sequence with label ---------------------------------------------
  proteinsequence <- reactive({sequences[sequences$Gene_ID == query(),]})
  uniprot_code    <- reactive({sequences[sequences$Gene_ID == query(),7]}) 
  labelposition   <- reactive({as.numeric(subset_two()$pos_aa)})
  labelposition_y <- reactive({rep(0, length(labelposition()))})
  label_content   <- reactive({paste0(subset_two()$ref_aa, subset_two()$pos_aa, subset_two()$alt_aa, collapse=",")})
  label_content2  <- reactive({str_split_fixed(label_content(), ',', Inf)})
  label_content3  <- reactive({as.data.frame(t(label_content2()))})
  
  # gene count --------------------------------------------------------------
  #############
  geneNAME0  <- reactive({subset(subset_two(), !grepl("subset of |covers", subset_two()$GeneSymbol, ignore.case=TRUE)) })
  geneNAME1  <- reactive({paste0(geneNAME0()$GeneSymbol, collapse = ";") })
  geneNAME   <- reactive({ gsub("-","Intergenic", geneNAME1())}) 
  #############
  gene0    <- reactive({as.data.frame(t(str_split(geneNAME(), pattern = ";", simplify = T)), stringsAsFactors = T) })
  gene1    <- reactive({as.data.frame(table(gene0()$V1)) })
  gene2    <- reactive({head( gene1()[order(- gene1()$Freq),], 10) })
  gene3    <- reactive({gene1()[order(- gene1()$Freq),] })
  
  
  output$displayngenes <- renderText({
    if (query() %in% genes$gene) {
      return(paste0("Protein Mapping"))
    } else {
      return(paste0(nrow(gene1()), " genes"))
    }
  })
  
  #only consequences which can be plotted
  subset_twoconsequence_A <- reactive({ gsub("3-UTR",NA, subset_two()$consequence)})
  subset_twoconsequence_B <- reactive({ gsub("5-UTR",NA, subset_twoconsequence_A())})
  subset_twoconsequence_C <- reactive({ gsub("Intronic",NA, subset_twoconsequence_B())})
  subset_twoconsequence_D <- reactive({ gsub("Mitocondrial",NA, subset_twoconsequence_C())})
  subset_twoconsequence_E <- reactive({ gsub("Non-coding",NA, subset_twoconsequence_D())})
  subset_twoconsequence_F <- reactive({ gsub("Splice-D/A",NA, subset_twoconsequence_E())})
  subset_twoconsequence_Fb <- reactive({ gsub("NA",NA, subset_twoconsequence_F())})
  subset_twoconsequence_G <- reactive({ gsub("Genomic",NA, subset_twoconsequence_Fb())})
  subset_twoconsequence_Gb <- reactive({ gsub("<NA>",NA, subset_twoconsequence_G())})  

  subset_twoconsequence_H <- reactive({ gsub("Frameshift","PTV", subset_twoconsequence_Gb())})
  subset_twoconsequence_I <- reactive({ gsub("In frame indel","PTV", subset_twoconsequence_H())})
  subset_twoconsequence_J <- reactive({ gsub("Stop gain","PTV", subset_twoconsequence_I())})
  subset_twoconsequence_K <- reactive({ as.factor(subset_twoconsequence_J())})
  
  
  output$GeneCountorProteinseq <- renderPlot({      
    if (query() %in% genes$gene) {      
      uniprot_code()  %>%
        drawProteins::get_features() %>%
        drawProteins::feature_to_dataframe() -> prot_data
      p <- ggplot(subset_two()) +
        geom_segment(mapping = aes(x= labelposition(),
                                   y= 1,
                                   xend = labelposition(),
                                   yend=3), show.legend = F) + 
        geom_point(aes(x= labelposition(), y=3,
                       #fill= subset_twoconsequence_K()),
                       colour= subset_twoconsequence_K(), show.legend = F),
                   size = 3, shape = 19) + 
        scale_colour_manual(breaks=c("Synonymous", "Missense", "PTV", NA),
                             values = c("Missense" = "grey34", "Synonymous" = "yellowgreen", "PTV" = "orangered1", NA))+
        # scale_fill_manual(breaks=c("Synonymous", "Missense", "PTV"),
        #                    values = c("Missense" = "grey", "Synonymous" = "green", "PTV" = "black", NA))+
        geom_point(aes(x=10, y= 4.75),size = 3, shape = 19, colour="orangered1") +
        geom_point(aes(x=10, y= 4.5),size = 3, shape = 19, colour="grey34") + 
        geom_point(aes(x=10, y= 4.25),size = 3, shape = 19, colour="yellowgreen") +
        
        geom_point(aes(x=10, y= 3.75), pch=21, size = 5, shape = 0) +
        geom_point(aes(x=10, y= 3.5), pch=21, size = 5, shape = 19, fill="yellow", colour="black") +

        geom_text(aes(x=20, y= 4.75), label ="Protein truncating variant (Stop-gain + Frameshift + Indel)", hjust = "left") +
        geom_text(aes(x=20, y= 4.5), label ="Missense", hjust = "left") +
        geom_text(aes(x=20, y= 4.25), label ="Synonymous", hjust = "left") +
        geom_text(aes(x=20, y= 3.75), label ="Uniprot domain", hjust = "left") +
        geom_text(aes(x=20, y= 3.5), label ="Phosphorlation site", hjust = "left") +
        
        # annotate("text", x= 10, y= 4.75, label = "Protein Truncating Variant") +
        # annotate("text", x= 10, y= 4.5, label = "Missense") +
        # annotate("text", x= 10, y= 4.25, label = "Synonimous") +
    
        xlim(1,proteinsequence()[[4]][1]) +
        theme_bw()+
        theme(panel.grid.minor=element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(),
              legend.title = element_blank(),
              axis.title = element_blank(),
              legend.position= "none") +
        ylim(0,5) +
        coord_cartesian(expand = F)
      p <- draw_chains(p, prot_data, label_chains = T)
      p <- draw_domains(p, prot_data)
      p <- draw_repeat(p, prot_data)
      p <- draw_motif(p, prot_data)
      p <- draw_phospho(p, prot_data, size = 6)      
      return(p)
    } else {      
      q <- ggplot(gene2(), mapping = aes(x = reorder(gene2()$Var1, gene2()$Freq), y = gene2()$Freq) ) +
        geom_bar(stat = "identity", fill = "#F04124", show.legend = F, width = 0.7) +
        ylab("Counts") +
        theme_light() +
        theme(title = element_blank(),
              axis.title.x=element_blank())
      return(q)
    }
  })
  
  # gene UI display decide --------------------------------------------------
  
  output$genefield <- renderUI({
    if (query() %in% genes$gene) {             
      fluidRow(column(width = 12, h6(paste0("Coding variant mapping and domain annotation for ",query())), align = "center",
                      withSpinner(plotOutput("GeneCountorProteinseq")))
      )
    } else {
      fluidRow(column(width = 8, h6("Top 10 genes associated"), align = "center", br(), 
                      withSpinner(plotOutput("GeneCountorProteinseq"))),
               column(width = 4, h6("Total genes associated"), align = "center",
                      wellPanel(DT::dataTableOutput("GeneCountTotal"),
                                style = "background-color: #ffffff;
                                border-color: #ffffff;
                                box-shadow: 0px 0px 0px #ffffff;
                                margin-bottom: 5px")
                      )
                      )
    }
    })
  
  # pheno count ------------------------------------------------------------
  phenoNAME <- reactive({ paste0(subset_two()$PhenotypeList, collapse = ";") })
  pheno0    <- reactive({ as.data.frame(t(str_split(phenoNAME(), pattern = ";", simplify = T)), stringsAsFactors = T) })
  pheno1    <- reactive({ as.data.frame(table(pheno0()$V1)) })
  pheno2    <- reactive({ head( pheno1()[order(- pheno1()$Freq),], 10) })
  pheno3    <- reactive({ pheno1()[order(- pheno1()$Freq),] })
  
  output$displayphenotypes <- renderText({
    return(paste0(nrow(pheno1()), " phenotypes"))
  })
  output$DiseaseCount <- renderPlot({
    ggplot(pheno2(), mapping = aes(x = reorder(stringr::str_wrap(pheno2()$Var1, 15), pheno2()$Freq), y = pheno2()$Freq) ) +
      geom_bar(stat = "identity", fill = "#E99002", show.legend = F, width = 0.7) +
      ylab("Counts") +
      theme_light() +
      theme(title = element_blank(),
            axis.title.x= element_blank())
  })
  
  # Table Final -----------------------------------------------------------------
  # Variantstable
  subset_two_table <- reactive({select(subset_two(), "GeneSymbol", "Type", "consequence", "ClinicalSignificance","review","PhenotypeList", "Name", "ref_aa", "alt_aa", "pos_aa") })
  subset_two_SNV <- reactive({subset_two_table()[subset_two_table()$Type == "SNV",]})
  subset_two_noSNV <- reactive({subset_two_table()[subset_two_table()$Type != "SNV",]})
  
  output$clinvartable1 <- DT::renderDataTable({
    DT::datatable(subset_two_SNV(), extensions = c('Scroller', 'Buttons'), 
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    'Table 1: ', htmltools::em('ClinVar single nucleotide variants (SNV) subset')),
                  rownames = F, colnames = c('Gene', 'Type', "Consequence", "Clinical significance","Review status", 'Phenotype', 'HGVS', "Ref_aa", "Alt_aa", "Pos_aa"), 
                  fillContainer = T,
                  options = list(dom = "t", 
                                 scrollY = 250,
                                 scroller = TRUE,
                                 pageLength = nrow(subset_two()),
                                 buttons = list(list(extend = 'collection',
                                                     buttons = c('csv', 'excel'),
                                                     text = 'Download'))
                                 ))
  })
  
  output$clinvartable2 <- DT::renderDataTable({
    DT::datatable(subset_two_noSNV(), extensions = c('Scroller', 'Buttons'),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    'Table 2: ', htmltools::em('ClinVar Copy number variants (CNV) subset')),
                  
                  rownames = F, colnames = c('Gene', 'Type', "Consequence", "Clinical significance","Review status",'Phenotype', 'HGVS', "Ref_aa", "Alt_aa", "Pos_aa"), 
                  fillContainer = T,
                  options = list(dom = "t", 
                                 scrollY = 250,
                                 scroller = TRUE,
                                 pageLength = nrow(subset_two()),
                                 buttons = list(list(extend = 'collection',
                                                     buttons = c('csv', 'excel'),
                                                     text = 'Download'))
                  ))
  })
  
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste0("SCV_", query(), "_" ,Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(subset_two_SNV(), file)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0("SCV_", query(), "_" ,Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(subset_two_noSNV(), file)
    }
  )
  
  # GeneCountTotaltable
  output$GeneCountTotal <- DT::renderDataTable({
    dat2 <- gene3()
    colnames(dat2) <- c("Gene", "Frequency")
    DT::datatable(dat2, 
                  #extensions = c('Scroller'), 
                  rownames = F, fillContainer = T,
                  options = list(dom = "t",
                                 scrollY = 350,
                                 scroller = TRUE,
                                 pageLength = nrow(gene1())))
  })
  # DiseaseTotaltable
  output$DiseaseCountTotal <- DT::renderDataTable({
    dat3 <- pheno3()
    colnames(dat3) <- c("Phenotype", "Frequency")
    DT::datatable(dat3, 
                  #extensions = c('Scroller'), 
                  rownames = F, 
                  fillContainer = T,
                  options = list(dom = "t",
                                 scrollY = 350,
                                 scroller = TRUE,
                                 pageLength = nrow(pheno1())
                  ))
  })
  # # Test
   #output$test <- renderPrint({
   #  print(nrow(clinicalgene_10()))
  # })
})
