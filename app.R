library(shiny)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Gene splicer"),
  
  # Sidebar with a slider input for number of bins 
  fluidRow(
    column(5,
           textAreaInput("cds_seq",
                      "Please input coding sequence (CDS)",
                     value = "  1 CCCTTTAAAG CCTTTAAAG"),
           textAreaInput("gene_seq",
                     "Please input the gene sequence (full genomic)",
                     value=" 1 GGCCCTTTAA AGGGCCCTTT AAAGGGCCCT TT"),
           textAreaInput("cds_reg",
                     "Please enter coding regions (copy and paste from TAIR, or enter numbers in sequence with a comma between each.",
                     value="coding_region	3-12	coding_region	 	16-24	 ")
           
    ),
    column(6,
           h3("The sequences flanking the splice point: "),
           tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
           verbatimTextOutput("seq"),
           h3("The position of insertions in the original CDS: "),
           tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
           verbatimTextOutput("ips"),
           h3("The spliced sequence with full stops inserted at splice points: "),
           tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
           verbatimTextOutput("outseq")
           #textOutput("search_seqs")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  library(tidyverse)
  
  reformat_seq <- function(string) {
    str_extract_all(string, "\\w") %>% 
      unlist() %>%
      str_extract_all("[A-Z]") %>%
      unlist() %>% str_flatten()
  }
  
  # reformatted_cds <- reactive({reformat_seq(input$cds_seq)})
  # #str_count(reformatted_cds) # correct
  # 
  # coding_regions <- reactive({input$cds_reg})
  # 
  # formatted_cdrs <- coding_regions %>%
  #   str_extract_all("([^0-9]|^)[0-9]+") %>%
  #   unlist() %>%
  #   str_extract_all("[0-9]+")
  # 
  # formatted_cdrs2 <- formatted_cdrs[2:(length(formatted_cdrs)-1)]
  # 
  # starts <- formatted_cdrs2[1:length(formatted_cdrs)%%2!=0] %>%
  #   map(as.numeric) %>%
  #   purrr::discard(function(x) identical(x, numeric(0)))
  # 
  # ends <- formatted_cdrs2[1:length(formatted_cdrs)%%2==0] %>%
  #   map(as.numeric) %>%
  #   purrr::discard(function(x) identical(x, numeric(0)))
  # 
  # gene_seq_f <- reactive({reformat_seq(input$gene_seq)}) 
  # 
  # prefixes <- list()
  # 
  # for(i in seq_along(starts)) {
  #   prefixes[[i]] <- str_sub(gene_seq_f, start=starts[[i]]-8, end=starts[[i]])
  # }
  # 
  # suffixes <- list()
  # 
  # for(i in seq_along(ends)) {
  #   suffixes[[i]] <- str_sub(gene_seq_f, start=ends[[i]], end=ends[[i]]+7)
  # }
  # 
  # search_sequences <- map2(prefixes, suffixes, str_c)
  # 
  # ips <- reformatted_cds %>% str_locate_all(unlist(prefixes))
  # 
  # ips <- map(ips, ~`[[`(.,2)) %>% unlist()
  # 
  # insert_dot <- function(string, insert_point) {
  #   output <- str_c(str_sub(string, start=1, end = insert_point), ".", str_sub(string, start=insert_point+1, end=str_count(string)))
  #   output
  # }
  # 
  # outseq <- reformatted_cds
  # for (i in seq_along(ips)) {
  #   outseq <- insert_dot(outseq, ips[[i]])
  # }
  
  output$seq <- renderPrint({
    reformatted_cds <- reformat_seq(input$cds_seq)
    #str_count(reformatted_cds) # correct
    
    coding_regions <- input$cds_reg 
    
    formatted_cdrs <- coding_regions %>%
      str_extract_all("([^0-9]|^)[0-9]+") %>%
      unlist() %>%
      str_extract_all("[0-9]+")
    
    formatted_cdrs2 <- formatted_cdrs[2:(length(formatted_cdrs)-1)]
    
    starts <- formatted_cdrs2[1:length(formatted_cdrs)%%2!=0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    ends <- formatted_cdrs2[1:length(formatted_cdrs)%%2==0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    gene_seq_f <- reformat_seq(input$gene_seq)
    
    prefixes <- list()
    
    for(i in seq_along(starts)) {
      prefixes[[i]] <- str_sub(gene_seq_f, start=starts[[i]]-8, end=starts[[i]])
    }
    
    suffixes <- list()
    
    for(i in seq_along(ends)) {
      suffixes[[i]] <- str_sub(gene_seq_f, start=ends[[i]], end=ends[[i]]+7)
    }
    
    prefixes <- map2(prefixes, "|", str_c)
    search_sequences <- map2(prefixes, suffixes, str_c)
    
    search_sequences
  })
  
  output$ips <- renderPrint({
    reformatted_cds <- reformat_seq(input$cds_seq)
    #str_count(reformatted_cds) # correct
    
    coding_regions <- input$cds_reg 
    
    formatted_cdrs <- coding_regions %>%
      str_extract_all("([^0-9]|^)[0-9]+") %>%
      unlist() %>%
      str_extract_all("[0-9]+")
    
    formatted_cdrs2 <- formatted_cdrs[2:(length(formatted_cdrs)-1)]
    
    starts <- formatted_cdrs2[1:length(formatted_cdrs)%%2!=0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    ends <- formatted_cdrs2[1:length(formatted_cdrs)%%2==0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    gene_seq_f <- reformat_seq(input$gene_seq)
    
    prefixes <- list()
    
    for(i in seq_along(starts)) {
      prefixes[[i]] <- str_sub(gene_seq_f, start=starts[[i]]-8, end=starts[[i]])
    }
    
    suffixes <- list()
    
    for(i in seq_along(ends)) {
      suffixes[[i]] <- str_sub(gene_seq_f, start=ends[[i]], end=ends[[i]]+7)
    }
    
    search_sequences <- map2(prefixes, suffixes, str_c)
    
    ips <- reformatted_cds %>% str_locate_all(unlist(prefixes))
    
    ips <- map(ips, ~`[[`(.,2)) %>% unlist()
    
    ips-1
  })
  
  output$outseq <- renderPrint({
    reformatted_cds <- reformat_seq(input$cds_seq)
    #str_count(reformatted_cds) # correct
    
    coding_regions <- input$cds_reg 
    
    formatted_cdrs <- coding_regions %>%
      str_extract_all("([^0-9]|^)[0-9]+") %>%
      unlist() %>%
      str_extract_all("[0-9]+")
    
    formatted_cdrs2 <- formatted_cdrs[2:(length(formatted_cdrs)-1)]
    
    starts <- formatted_cdrs2[1:length(formatted_cdrs)%%2!=0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    ends <- formatted_cdrs2[1:length(formatted_cdrs)%%2==0] %>%
      map(as.numeric) %>%
      purrr::discard(function(x) identical(x, numeric(0)))
    
    gene_seq_f <- reformat_seq(input$gene_seq)
    
    prefixes <- list()
    
    for(i in seq_along(starts)) {
      prefixes[[i]] <- str_sub(gene_seq_f, start=starts[[i]]-8, end=starts[[i]])
    }
    
    suffixes <- list()
    
    for(i in seq_along(ends)) {
      suffixes[[i]] <- str_sub(gene_seq_f, start=ends[[i]], end=ends[[i]]+7)
    }
    
    search_sequences <- map2(prefixes, suffixes, str_c)
    
    ips <- reformatted_cds %>% str_locate_all(unlist(prefixes))
    
    ips <- map(ips, ~`[[`(.,2)) %>% unlist()
    
    insert_dot <- function(string, insert_point) {
      output <- str_c(str_sub(string, start=1, end = insert_point), ".", str_sub(string, start=insert_point+1, end=str_count(string)))
      output
    }
    
    outseq <- reformatted_cds
    
    ips_new <- ips + 0:(length(ips)-1) # each time you insert a dot, you need to increase the IP by one.
    
    for (i in seq_along(ips_new)) {
      outseq <- insert_dot(outseq, ips_new[[i]])
    }
    
    outseq
    })
  
  #output$search_seqs <- renderText({search_sequences})
  
}

# Run the application 
shinyApp(ui = ui, server = server)
