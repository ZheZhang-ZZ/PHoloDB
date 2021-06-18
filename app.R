library(shinydashboard)
library(shiny)
library(plotly)
library(tidyverse)
library(ggpubr)
library(shinyjs)
library(shinyWidgets)
library(data.table)
library(waiter)
library(ggsci)
library(leaflet)
library(leaflet.extras) 
#library(shinythemes)
#library(shinydashboardPlus)

# directory
gen_dir <- "example_data/geno_data"
phe_dir <- "example_data/phe_data"
meta_dir <- "example_data/meta_data"

# read data
farms<-c("SHD","GXL","GXH")
farms_title<-list(
  SHD="SHD",
  `GXL (low IMF line)`="GXL",
  `GXH (high IMF line)`="GXH"
)
farms_name<-c("SHD","GXL (low IMF line)","GXH (high IMF line)")

meta_farms<-c("GXL","GXH")
meta_farms_title<-list(
  `GXL (low IMF line)`="GXL",
  `GXH (high IMF line)`="GXL"
)

meta_farms_name<-c("GXL (low IMF line)","GXH (high IMF line)")

traits<-c("Birth_weight","Day_to_115kg", "Backfat_thickness" ,"Loin_muscle_area", "Lean_percentage")
traits_name<-c("Birth weight", "Days to 115 kg", "Backfat thickness",
               "Loin muscle area", "Lean percentage")
traits_title<-list(Birth_weight="Birth weight (kg)",
                   Day_to_115kg="Days to 115 kg (d)",
                   Backfat_thickness="Backfat thickness (mm)",
                   Loin_muscle_area="Loin muscle area (cm<sup>2</sup>)",
                   Lean_percentage="Lean percentage (%)")

chrs<-paste0("Chr",18)
gen_size<-map_chr(chrs,function(x){
  return(paste0(x,":         ",round(file.info(paste0(gen_dir,"/",x,".wgs.snp.zip"))$size/1024^2,1), "Mb"))
})

## phenotype
phe_files<-c(paste0(phe_dir,"/SHD_phe.csv"),
             paste0(phe_dir,"/GXL_phe.csv"),
             paste0(phe_dir,"/GXH_phe.csv"))

phe<-map(phe_files,fread)
names(phe)<-farms

## genotype
### sample size
all_geno_info<-read.table(paste0(gen_dir, "/all_ind_geno_info.txt"),
                          h=T,sep="\t")
all_geno_info$Farm<-factor(all_geno_info$Farm,levels=farms)
all_geno_info$Platform<-factor(all_geno_info$Platform,levels=c('LC','ULC','CHIP'))

### PCA
pca_dat<-read.table(paste0(gen_dir, "/pca_plot.txt"),h=T,sep="\t")
pca_dat$GRP<-factor(pca_dat$GRP,levels=farms)

### sequence summary
gen_seq_sum_tab<-read.table(paste0(gen_dir,"/seq_ind_geno_info.txt"),
                            h=T,sep="\t")

### SNP distribution file
dist_snp<-read.table(paste0(gen_dir,"/snp_dist.txt"),h=T,sep="\t")
dist_snp$Chromosome<-paste0("Chr",dist_snp$Chromosome)
dist_snp$Chromosome<-factor(dist_snp$Chromosome,levels=chrs)

## Metagenomic data

### sequencing summary
meta_seq_sum_tab<-read.table(paste0(meta_dir,"/summaryBase_meta.txt"),
                            h=T,sep="\t")
colnames(meta_seq_sum_tab)<-c("ID", "Breed","Line","RawReads", "CleanReads",
              "Q20 (%)", "Q30 (%)")

## pylum redundance
phylum_red<-read.table(paste0(meta_dir,"/phylum_red.txt"),h=T,sep="\t")
phylum_legend<-phylum_red %>%
  group_by(micro) %>%
  summarise(SUM=sum(redundance)) %>%
  arrange(desc(SUM)) %>%
  select(micro)
phylum_legend<-phylum_legend[['micro']]
phylum_legend<-as.character(phylum_legend)

phylum_legend<-phylum_legend[!phylum_legend=='Others']
phylum_legend<-append(phylum_legend,'Others')
phylum_red$micro<-factor(phylum_red$micro,levels=phylum_legend)
phylum_red$grp<-factor(phylum_red$grp,levels=meta_farms)

## metagenomic pca data
meta_pca_dat<-read.table(paste0(meta_dir, "/pca_plot.txt"),h=T,sep="\t")
meta_pca_dat$Line<-factor(meta_pca_dat$Line,levels=meta_farms)

### Map info
#load("data/mapdata.rData") # China map info with Taiwan
map_info<-read.table("example_data/map_info.txt",h=T,sep="\t") # population map info
map_info$Lat <- jitter(map_info$Lat, factor = 0.0001)
map_info$Lng <- jitter(map_info$Lng, factor = 0.0001)
rgn<-c(73.60226, 15.77538, 134.77258, 53.56944) # China view region
carto = "http://a.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png"

# Home panel
homePanel<-tabPanel(
  title = "Home",
  
  # Title band
  fluidRow(
    width=12,
    style = "margin-top:-2em",
    
    box(title = "", 
        height = 230, 
        background="blue",
        width="100%",
        htmlOutput("text"),
        tags$head(
          tags$style("#text{color: white;
                            text-align:center;
                            font-size: 50px;
                            font-style: bold;
                            position:relative;
                            top:50px;
                            font-family: Sans-serif;
                            }"
          )
        )
    )
  ),
  
  # inline white space
  div(class = "inlay", style = "height:15px;width:100%;background-color:ghostwhite;"),
  
  # Map and info boxes
  fluidRow(
    column(
      width = 11, offset = 1,
      fluidRow(
        column( width = 8,
                
                # Sasmple distribution map ------------------------------------------------------
                box(title = "Samples distribution",
                    status = "primary",
                    solidHeader = TRUE,
                    width = NULL,
                    leafletOutput("map", height = 410))
        ),
        
        column( width = 3,
                tags$head(tags$style(HTML(".small-box {width:110%}"))),
                # valueBox1 ----------------------------------------------------------------
                fluidRow(
                  valueBox("3 different lines",
                           subtitle="Within 1 pig breed",
                           icon=icon("map-marker"),
                           color='red',width=12)
                  ),
                
                # valueBox2 ----------------------------------------------------------------
                fluidRow(
                  valueBox(
                    paste0(ncol(phe[[1]])-6, " traits"),
                    subtitle = format(paste0(18931," pigs with phenotype"),
                                      big.mark = " ",
                                      scientific = FALSE),
                    icon = icon("chart-pie"),
                    color = "yellow",
                    width = 12
                  )),
                
                # valueBox3 ----------------------------------------------------------------
                fluidRow(
                  # genotype info box
                    valueBox(
                      paste0(length(unique(all_geno_info$ID)), " genomes"),
                      subtitle = format(paste0(sum(dist_snp$Count)," high-qulity SNPs are available"),
                                        big.mark = " ",
                                        scientific = FALSE),
                      icon = icon("dna"),
                      color = "olive",
                      width = 12
                    )),
  
                # valueBox4 ----------------------------------------------------------------
                fluidRow(
                  valueBox(
                    paste0(47, " metagenomes"),
                    subtitle = "Microbiome down to species level is available",
                    icon = icon("virus"),
                    color = "maroon",
                    width = 12
                  ))
        ))
      )
  )
)

# Visualization panel
visualPanel<-tabPanel(
  title = "Data visualization",
  
  # Phenotype data
  fluidRow(
    column(
      width = 11, offset = 1,
      
      fluidRow(
        column(
          width = 12,
          div(class = "inlay", style = "height:100px;width:100%;background-color: ghostwhite;"),
          div(style="display:inline-block",
              actionButton(
                inputId = "pheBttn",
                label = "Phenotypic Data",
                #width = "100%",
                class = "btn-info",
                style='font-size:150%;width:210px'
                #style="color: #fff; background-color: #337ab7; border-color: light-blue"
              )
          ),
          div(class = "inlay", style = "height:15px;width:100%;background-color: ghostwhite;")
        )
      ),
      
      fluidRow(
        div(style="width:1430px;",
          box(
            width=2,
            title = "Settings",
            status = "primary",
            solidHeader = TRUE,
            #collapsible = TRUE,
            height=550,
            br(),
            awesomeRadio(
              inputId = "population",
              label = "Choose line:", 
              choices = c(farms_name,'ALL'),
              selected = "SHD",
              status = "warning"
            ),
            br(),
            awesomeRadio(
              inputId = "traits",
              label = "Choose trait:", 
              choices = traits_name,
              selected = "Backfat thickness",
              status = "danger"
            ),
            br(),
            awesomeRadio(
              inputId = "fix",
              label = "Choose factors:", 
              choices = c("Parity","Sex", "Herd"),
              selected = "Parity",
              status = "success"
            ),
            br()
          )
        ),
        uiOutput("pheBoxes")
      )
    )
  ),
  
  # Genotype data
  fluidRow(
    column(
      width = 11, offset = 1,
      fluidRow(
        column(
          width = 12,
          div(style="display:inline-block",
              actionButton(
                inputId = "pheBttn",
                label = "Genotypic Data",
                width = 210,
                class = "btn-success",
                style='font-size:150%'
                #style="color: #fff; background-color: #337ab7; border-color: light-blue"
              )
          ),
          div(class = "inlay", style = "height:15px;width:100%;background-color: ghostwhite;")
        ) 
      ),
      
      fluidRow(
        div(style="width:1430px;",
          box(
            width=2,
            title = "Settings",
            status = "primary",
            solidHeader = TRUE,
            #collapsible = TRUE,
            height=500,
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            awesomeRadio(
              inputId = "genItem",
              label = "Item to be shown:", 
              choices = c("Sequencing summary",
                          "Sample size",
                          "SNP distribution",
                          "PCA result"),
              selected = "Sequencing summary",
              status = "success"
            )
          )
        ),
        uiOutput('genBoxes')
      )
    )
  ),
  
  # Microbial data
  fluidRow(
    column(
      width = 11, offset = 1,
      fluidRow(
        column(
          width = 12,
          div(style="display:inline-block",
              actionButton(
                inputId = "microBttn",
                label = "Microbiome Data",
                width = 210,
                class = "btn-danger",
                style='font-size:150%'
                #style="color: #fff; background-color: #337ab7; border-color: light-blue"
              )
          ),
          div(class = "inlay", style = "height:15px;width:100%;background-color: ghostwhite;")
        ) 
      ),
      
      fluidRow(
        div(style="width:1430px;",
            box(
              width=2,
              title = "Settings",
              status = "primary",
              solidHeader = TRUE,
              #collapsible = TRUE,
              height=500,
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              awesomeRadio(
                inputId = "microItem",
                label = "Item to be shown:", 
                choices = c("Sequencing summary",
                            "Sample size",
                            "Taxonomic profile",
                            "PCA result"),
                selected = "Taxonomic profile",
                status = "danger"
              )
            )
        ),
        uiOutput('microBoxes')
      )
    )
  ),
  div(class = "inlay", style = "height:30px;width:100%;background-color: ghostwhite;")
)

# Download panel
downloadPanel<-tabPanel(
  title = "Data download",
  
  div(class = "inlay", style = "height:100px;width:100%;background-color: ghostwhite;"),
  
  fluidRow(
    column(
      width = 11, offset = 1,
      
      fluidRow(
        
        # phenotypic data download ---------------------------------------------------
        column(
          width = 4,
         div(style="width:150%;",
             box(
               title = "Phenotypic data",
               status = "primary",
               solidHeader = TRUE,
               height=560,
               br(),
               awesomeRadio(
                 inputId = "pheDownloadPop",
                 label = "Choose line:", 
                 choices = farms_name,
                 selected = "SHD",
                 status = "primary"
               ),
               br(),
               br(),
               br(),
               br(),
               br(),
               awesomeRadio(
                 inputId = "pheDownloadType",
                 label = "Choose data type:", 
                 choices = c("Phenotype", "Pedigree"),
                 selected = "Phenotype",
                 status = "danger"
               ),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
                downloadBttn(outputId="download_phe",
                            style="jelly",
                            color="success")
             )
          )
        ),
        
        # genotypic data download ---------------------------------------------------
        column(
          width = 4,
          div(style="width:150%;",
              box(
                title = "Genotypic data",
                status = "success",
                solidHeader = TRUE,
                height=560,
                br(),
                awesomeRadio(
                  inputId = "genDownloadType",
                  label = "Choose data type:", 
                  choices = c("Imputed WGS SNPs", "Chip SNPs"),
                  selected = "Imputed WGS SNPs",
                  status = "danger"
                ),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                selectInput("download_gen_chr", "Choose chromosome:",
                            choices = gen_size),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                downloadBttn(outputId="download_gen",
                             style="jelly",
                             color="success"),
                br()
              )
          )
        ),
        
        # microbial data
        column(
          width = 4,
          div(style="width:150%;",
            box(
              #width=2,
              title = "Microbial data",
              status = "danger",
              solidHeader = TRUE,
              height=560,
              br(),
              awesomeRadio(
                inputId = "microDownloadPop",
                label = "Choose line:", 
                choices = meta_farms_name,
                selected = meta_farms_name[1],
                status = "primary"
              ),
              br(),
              awesomeRadio(
                inputId = "microDownloadType",
                label = "Choose data type:", 
                choices = c("Relative abundance",
                            "Count"),
                selected = "Relative abundance",
                status = "warning"
              ),
              br(),
              awesomeRadio(
                inputId = "microDownloadClass",
                label = "Choose level:", 
                choices = c("Phylum","Class","Order","Family",
                            "Genus","Species"),
                selected = "Phylum",
                status = "danger"
              ),
              br(),
              downloadBttn(outputId="download_micro",
                           style="jelly",
                           color="success")
            )
          )
        )
        
      )
    )
  ),
  div(class = "inlay", style = "height:30px;width:100%;background-color: ghostwhite;")
)

document_panel<-tabPanel(
  title = "About PHoloDB",
  
  div(class = "inlay", style = "height:100px;width:100%;background-color: ghostwhite;"),
  
  # doc <- tags$html(
  #   tags$body(
  #     p(strong('P',style='color:#bcbd22',.noWS="after"),'ig',
  #       strong('H',style='color:#ff7f0e',.noWS="after"),'olog',
  #       strong('eno',style='color:#8c564b',.noWS="outside"), 'mic',
  #       strong('D',style='color:#d62728',.noWS="after"),'ata',
  #       strong('B',style='color:#9467bd',.noWS="after"),'ase',
  #       style='font-size:40px;
  #              position:relative;
  #              top:20px;
  #              font-family: Sans-serif;'
  #     )
  #   )
  # )
  
  fluidRow(
    column(
      width = 8, offset = 2,
      box(
        title = NULL, status = "primary", width = NULL,
        
        tags$div(
          class = "document",
          div(tags$b("Introduction", style='font-size:20px;')),
          br(),
          tags$div("Welcome to Pig Hologenomic Data Base (PHoloDB) v1.0 !"),
          br(),
          tags$div(tags$a("Hologenome", href = "https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0457-9"),
                   " means the sum of host genome and associated microbiome.
                    It is a revolutionary concept proposed over ten years ago, and has provided a novel viewpoint
                    to understand the formation, development, evolution, and characteristics of distinct biological entity. 
                    Pig is one of the most important meat suppliers and model animals for human beings. Therefore,
                    it will be of great importance to reveal the genetic mechanisms of pigs' phenotypes from hologenomic perspective."),
          br(),   
          tags$div("In order to help global scientists to make use of more integrated data to study the chracteristics of pigs using
                    hologenomic methods and validate novel methods that can be used in pig breeding or genetic mechanism 
                    dissection, we develop PHoloDB to curate data, inluding but not limited to pig and its carried microorganism's genome,
                    transcriptome, epigenome, phenotype, microbiome, metabolome, etc."),
          br(),         
          tags$div("In the present version (1.0) of PHoloDB, we include the phenotype, genotype and microbiome data from three different pig
                    lines. In the next version, we will further enlarge the sample size and also add the uploading function to make PHoloDB
                    become a database that help global scientists to curate the pig hologenomics associated data.
                   "),
          br(),
          
          div(tags$b("Function", style='font-size:20px;')),
          br(),
          tags$div("Based on the curated data of version 1.0, PHoloDB at present include the following functions:"
          ),
          br(),
          tags$ul(
            tags$li("Visualization of different omics data"),
            tags$li("Data download and the explanation of the data structure can be found ", 
                    tags$a("here", href = "https://github.com/ZheZhang-ZZ/PHoloDB")
            )
          ),
          br(),
          tags$div("In the next version, we will include more data and add data uploading functions."
          ),
          br(),
          
          div(tags$b("Framework and source code", style='font-size:20px;')),
          br(),
          tags$div("PHoloDB is nearly 100% developed by ",
                   tags$a("Shiny", href = "https://cran.r-project.org/web/packages/shiny/index.html"), 
                   " R package and its friend packages. The source code of PHoloDB can be found ",
                   tags$a("here", href = "https://github.com/ZheZhang-ZZ/PHoloDB"),"."
          ),
          br(),
          
          div(tags$b("Citation", style='font-size:20px;')),
          br(),
          tags$div("If you use any of the data or function of PHoloDB, please cite at least one of the
                   following publications:"),
          br(),
          tags$ul(
            "Comming soon..."
          ),
          br(),
          
          div(tags$b("Contact info", style='font-size:20px;')),
          br(),
          tags$div("PHoloDB was conceived by Yuchun Pan (panyuchun1963@aliyun.com) and
                   developed by Zhe Zhang (zhe_zhang@zju.edu.cn) ",
                   tags$a(class = "btn btn-default", icon("github"), "@Zhe", 
                          href = "https://github.com/ZheZhang-ZZ", 
                          style = "background-color: #24292E; color:#FFFFFF;"),
                   "and Zhenyang Zhang (zhangzy1995@aliyun.com)",
                   tags$a(class = "btn btn-default", icon("github"), "@Zhenyang", 
                          href = "https://github.com/zhangzhenyang", 
                          style = "background-color: #24292E; color:#FFFFFF;")
          ),
          br(),
          tags$div("You can also try our other online databases and tools in this website: ", 
                   tags$a("alphaindex.zju.edu.cn", href = "http://alphaindex.zju.edu.cn/"),
                   "."),
          br()
        )
        
      )
    )
    
  ),
  
  div(class = "inlay", style = "height:30px;width:100%;background-color: ghostwhite;")
)


ui<-navbarPage(
  title="Pig Hologenomic Data Base (PHoloDB)",
  
  position = 'fixed-top',
  header = tagList(
    useShinydashboard(),
    useShinydashboardPlus(),
    use_waiter(),
    useShinyjs(),
    setBackgroundColor(color = c("ghostwhite"))
  ),
  theme = shinythemes::shinytheme(theme = "flatly"),
  
  
  homePanel,
  visualPanel,
  downloadPanel,
  document_panel,
  
  footer = div(class="footer",
               includeHTML("footer.html")
  )
  
)

server<-function(input,output,session){
  
  # Download Panel -----------------------------------------------------------------------
  
  ## phenotypic data downlaod
  phe_to_download<-reactive({
    if(input$pheDownloadType=="Pedigree"){
      dat<-paste0(phe_dir,"/",
                  farms_title[[input$pheDownloadPop]],"_ped.csv")
    }else{
      dat<-paste0(phe_dir,"/",
                  farms_title[[input$pheDownloadPop]],"_phe.csv")
    }
    return(dat)
  })
  
  output$download_phe<-downloadHandler(
    filename = function() {
      if(input$pheDownloadType=="Phenotype"){
        return(paste0(farms_title[[input$pheDownloadPop]], "_phenotype.csv"))
      }else{
        return(paste0(farms_title[[input$pheDownloadPop]], "_pedigree.csv"))
      }
    },
    content = function(file) {
      file.copy(phe_to_download(), file)
    })
  
  ## Genotypic data download
  gen_file <- reactive({
    if (input$genDownloadType=='Chip SNPs') {
      return(paste0(gen_dir,"/chip.raw.snp.zip"))
    }else{
      return(paste0(gen_dir,"/",chrs[which(gen_size==input$download_gen_chr)],".wgs.snp.zip"))
    }
  })

  observeEvent(input$genDownloadType,
               {
                 if (input$genDownloadType!='Chip SNPs') {
                   shinyjs::enable("download_gen_chr")
                 }else{
                   shinyjs::disable("download_gen_chr")
                 }
               })

  output$download_gen<-downloadHandler(
    filename = function() {
      file_name<-unlist(str_split(gen_file(),'\\/'))
      return(file_name[length(file_name)])
    },
    content = function(file) {
      file.copy(gen_file(), file)
    },
    contentType = "application/zip"
  )
  
  ## Microbial data download
  micro_file <- reactive({
    if (input$microDownloadType=='Relative abundance') {
      return(paste0(meta_dir,"/",meta_farms_title[[input$microDownloadPop]],"_abundance_",input$microDownloadClass,".txt"))
    }else{
      return(paste0(meta_dir,"/",meta_farms_title[[input$microDownloadPop]],"_count_",input$microDownloadClass,".txt"))
    }
  })
  
  output$download_micro<-downloadHandler(
    filename = function() {
      file_name<-unlist(str_split(micro_file(),'\\/'))
      return(file_name[length(file_name)])
    },
    content = function(file) {
      file.copy(micro_file(), file)
    }
  )
  
  # Home tabPanel ------------------------------------------------------------------------

  ## title band
  doc <- tags$html(
    tags$body(
      p(strong('P',style='color:#bcbd22',.noWS="after"),'ig',
        strong('Hol',style='color:#ff7f0e',.noWS="after"),'ogen',
        strong('o',style='color:#8c564b',.noWS="outside"), 'mic',
        strong('D',style='color:#d62728',.noWS="after"),'ata',
        strong('B',style='color:#9467bd',.noWS="after"),'ase',
        style='font-size:40px;
               position:relative;
               top:20px;
               font-family: Sans-serif;'
      )
    )
  )
  
  output$text <- renderUI(
    doc
  )
  
  ## Sample distribution map
  
  output$map<-renderLeaflet({
    map_info %>%
      leaflet() %>% addTiles(carto) %>% clearShapes() %>%
      #addProviderTiles("CartoDB", group = "Carto") %>%
      flyToBounds(rgn[1], rgn[2], rgn[3], rgn[4])%>%
      #addPolygons(data = mapdata, fill = F, color = 'red') %>%
      setView(lng = 112.939003, lat = 28.228001, zoom = 4) %>%
      addMarkers(
        lng=~Lng,
        lat=~Lat,
        popup=~paste0("Breed: ", Breed, "<br/>", 
                      "Line: ", Line, "<br/>",
                      "Sample size: ", format(Size, big.mark=" ",scientific=FALSE))) %>%
        addResetMapButton() 
  })
  
  # Microbiomic data tabPanel  ----------------------------------------------------------------
  
  output$microBoxes<-renderUI({
    switch(input$microItem,
           `Sequencing summary` = box(
             title = "Summary of metagenomic sequencing data", status = "primary",
             collapsible = FALSE,
             DT::dataTableOutput("meta_seq_sum_tab",width="100%"),
             width = 9,height = 500
           ),
           `Sample size` = box(
             title = "Sample size of metagenomic sequenced pigs", status = "primary",
             collapsible = FALSE,
             plotlyOutput("meta_sample_size_plot", height = 430),
             width = 9, height = 500
           ),
           `Taxonomic profile` = tabBox(title = 'Taxonomic classification profiling',
                                       width=9,id = "tabsetTax", height = 500,
                                       tabPanel("Abundance plot", plotlyOutput("abc_plot", height = 430)),
                                       tabPanel("Cladogram",imageOutput("cladegram"))),
           `PCA result` = box(
             title = "PCA of metagenomic sequenced pigs by Bray-Curtis distance", status = "primary",
             collapsible = FALSE,
             plotlyOutput("meta_pca_plot", height = 430),
             width = 9, height = 500
           )
    )
  })
  
  ## Summary table output
  
  output$meta_seq_sum_tab <- DT::renderDataTable(meta_seq_sum_tab, 
                                                options = list(pageLength = 7,
                                                               lengthChange = FALSE))
  
  ## Sample size plot
  meta_sample_size_plot<-reactive({
    req(input$microItem)
    waiter<-Waiter$new(id='meta_sample_size_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    p<-meta_seq_sum_tab %>%
      group_by(Line) %>%
      summarise(Size=n()) %>%
      plot_ly() %>%
      add_bars(x = ~Line, 
               y = ~Size,
               color = ~Line, 
               colors="#d62728",
               opacity=0.9) %>%
      layout(
        yaxis = list(
          title = "Sample size"
        ),showlegend = FALSE)
    return(p)
  })
  
  output$meta_sample_size_plot<-renderPlotly(meta_sample_size_plot())
  
  ## Phylum abundance plot
  abc_plot<-reactive({
    req(input$microItem)
    waiter<-Waiter$new(id='meta_sample_size_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    
    red_plot<-ggplot(phylum_red,aes(x = ID, y=redundance, fill = micro)) + 
      geom_bar(stat='identity',width=0.7,
               position = position_fill(reverse = TRUE))+
      facet_grid(. ~ grp, scales = "free", space='free')+
      theme_pubr()+
      scale_y_continuous(breaks=seq(from=0,to=1,by=0.1),
                         labels=seq(from=0,to=100,by=10))+
      scale_fill_manual(values=pal_d3("category20c")(15))+
      theme(axis.ticks.x = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            text = element_text(size=20),
            axis.text.x = element_blank())+ylab("Relative abundance (%)")
    
    fig <- ggplotly(red_plot) %>%
      layout(
        legend = list(
          orientation = "h",   # show entries horizontally
          xanchor = "center",  # use center of legend as anchor
          x = 0.5,y=-0.105),
        font=list(
          family = "sans serif",
          size = 15
        ))
    return(fig)
  })
  
  output$abc_plot<-renderPlotly(abc_plot())
  
  ## Cladegram
  output$cladegram<-renderImage({
    list(
      src = paste0(meta_dir,"/graphlan.jpg"),
      contentType = "image/jpeg",
      width = 600,
      height = 400,
      style="display: block; margin-left: auto; margin-right: auto;"
    )
  }, deleteFile = FALSE)
  
  ## PCA plot
  meta_pca_plot<-reactive({
    req(input$microItem)
    waiter<-Waiter$new(id='meta_pca_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    
    p<-ggplot(meta_pca_dat,aes(x=PC1,y=PC2,colour=Line,shape=Line,text=ID))+
      geom_point(size=3,alpha=0.8)  +
      scale_colour_manual(values = c("#00AFBB", "#FC4E07"))+
      scale_shape_manual(values=c(15,19)) +
      #scale_x_continuous(breaks=seq(from=0,to=0.08,by=0.02))+
      theme_bw()+
      theme(legend.title=element_blank(),
            legend.position = c(0.85,0.8),
            text = element_text(size=20),
            axis.title.x = element_text(vjust=-1.1),
            axis.text = element_text(size=20,color='black'),
            legend.background = element_rect(fill=alpha('white', 0.0))) +
      xlab('PC1 (22.2.00%)')+ylab('PC2 (13.4%)')
    
    fig <- ggplotly(p,tooltip='ID') %>%
      layout(legend = list(
        x=0.5,y=-0.4,
        orientation = "h",   # show entries horizontally
        xanchor = "center" # use center of legend as anchor
      ),
      font=list(
        family = "sans serif",
        size = 20
      ))
    return(fig)
  })
  
  output$meta_pca_plot<-renderPlotly(meta_pca_plot())
  
  # Genotype data tabPanel  ----------------------------------------------------------------
  
  output$genBoxes<-renderUI({
    switch(input$genItem,
           `Sequencing summary` = box(
             title = "Summary of sequencing data", status = "primary",
             collapsible = FALSE,
             DT::dataTableOutput("gen_seq_sum_tab",width="100%"),
             width = 9,height = 500
           ),
           `Sample size` = box(
             title = "Sample size of different genotyping platforms", status = "primary",
             collapsible = FALSE,
             plotlyOutput("gen_sample_size_plot", height = 430),
             width = 9, height = 500
           ),
           `SNP distribution` = tabBox(title = 'Imputed WGS SNP density and distribution',
                                       width=9,id = "tabsetSNP", height = 500,
                                       tabPanel("Distribution plot", plotlyOutput("snp_dis_plot", height = 430)),
                                       tabPanel("Density plot",imageOutput("snp_dens_plot"))),
           `PCA result` = box(
             title = "PCA of all populations based on imputed SNP panel", status = "primary",
             collapsible = FALSE,
             plotlyOutput("gen_pca_plot", height = 430),
             width = 9, height = 500
           )
    )
  })
  
  ## sample_size_plot
  gen_sample_size_plot<-reactive({
    req(input$genItem)
    waiter<-Waiter$new(id='gen_sample_size_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    p<-all_geno_info %>%
      group_by(Farm,Platform) %>%
      summarise(Size=n()) %>%
      plot_ly() %>%
      add_bars(x = ~Farm, 
               y = ~Size,
               color = ~Platform, 
               colors = "Dark2",
               opacity=0.9) %>%
      layout(
        yaxis = list(
          title = "Sample size"
        ),
        legend = list(x = 100, y = 0.5))
    return(p)
  })
  
  output$gen_sample_size_plot<-renderPlotly(gen_sample_size_plot())
  
  ## PCA plot
  gen_pca_plot<-reactive({
    req(input$genItem)
    waiter<-Waiter$new(id='gen_pca_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    
    p<-ggplot(pca_dat,aes(x=PC1,y=PC2,colour=GRP,shape=GRP,text=ID))+
      geom_point(size=3,alpha=0.8)  +
      scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07"))+
      scale_shape_manual(values=c(17,15,19)) +
      scale_x_continuous(breaks=seq(from=0,to=0.08,by=0.02))+
      theme_bw()+
      theme(legend.title=element_blank(),
            legend.position = c(0.1,0.2),
            text = element_text(size=20),
            axis.title.x = element_text(vjust=-1.1),
            axis.text = element_text(size=20,color='black'),
            legend.background = element_rect(fill=alpha('white', 0.0),
                                             colour='black')) +xlab('PC1 (8.00%)')+ylab('PC2 (2.54%)')
    fig <- ggplotly(p,tooltip='ID') %>%
      layout(legend = list(
        x=0.8,y=0.1
      ),
      font=list(
        family = "sans serif",
        size = 20
      ))
    return(fig)
  })
  
  output$gen_pca_plot<-renderPlotly(gen_pca_plot())
  
  ## SNP distribution plot
  snp_dis_plot<-reactive({
    req(input$genItem)
    waiter<-Waiter$new(id='snp_dis_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    dist_snp %>%
      plot_ly() %>% 
      add_bars(x = ~Chromosome, y = ~Count) %>%
      layout(
        xaxis = list(
          title = "Chromosome"
        ),
        yaxis = list(
          title = "Number of SNPs"
        ),
        
        font=list(
          family = "sans serif",
          size = 20
        ))
  })
  output$snp_dis_plot<-renderPlotly(snp_dis_plot())
  
  ## SNP density plot
  output$snp_dens_plot<-renderImage({
    list(
      src = paste0(gen_dir,"/SNP-Density.Pval.jpg"),
      contentType = "image/jpeg",
      width = 800,
      height = 440,
      style="display: block; margin-left: auto; margin-right: auto;"
    )
  }, deleteFile = FALSE)
  
  ## Summary table output
  
  output$gen_seq_sum_tab <- DT::renderDataTable(gen_seq_sum_tab, 
                                            options = list(pageLength = 7,
                                                           lengthChange = FALSE))
  
  
  # Phenotype data tabPanel  ----------------------------------------------------------------
  observeEvent(input$population,{
    if(input$population=='ALL'){
      shinyjs::disable('fix')
    }else{
      shinyjs::enable('fix')
    }
  })
  
  output$pheBoxes<-renderUI({
    if(input$population!='ALL'){
      tabBox(
        title = paste0("Plots of traits for ", toupper(farms_title[[input$population]]), " population"),
        width=9,
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1", height = 550,
        tabPanel("Distribution", 
                 plotlyOutput("trait_plot", height = 400)
        ),
        tabPanel("Sample size",
                 plotlyOutput("size_plot", height = 400)
        )
      )
    }else{
      tabBox(
        title = "Plots of traits for all populations",
        width=9,
        id = "tabset1", height = 550,
        tabPanel("Distribution", 
                 plotlyOutput("trait_plot", height = 400)
        ),
        tabPanel("Sample size",
                 plotlyOutput("size_plot", height = 400)
        )
      )
    }
  })
  
  phe_farm<-reactive({
    if(input$population!='ALL'){
      farm<-farms_title[[input$population]]
      return(phe[[toupper(farm)]])
    }else{
      phe2<-rbindlist(phe,idcol = "Farm")
      phe2$Farm<-factor(phe2$Farm,levels=farms)
      return(phe2)
    }
  })
  
  dens_plots<-reactive({
    req(input$traits)
    waiter<-Waiter$new(id='trait_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    selected_trait<-traits[traits_name%in%input$traits]
    if(input$population!='ALL'){
      req(input$fix)
      map2(selected_trait, input$fix,function(x,y){
        p<-phe_farm() %>%
          plot_ly(
            x= .[[y]],
            y= .[[x]],
            split = .[[y]],
            type = 'violin',
            box=list(visible=T),
            meanline=list(visible=T)
          ) %>%
          layout(
            xaxis = list(
              title = y
            ),
            yaxis = list(
              title = traits_title[[x]]
              #zeroline = F
            ),
            font=list(
              family = "sans serif",
              size = 20
            ),
            showlegend=FALSE
          )
        return(p)
      })
    }else{
      map(selected_trait, function(x){
        p<-phe_farm() %>%
          plot_ly(
            x= .[['Farm']],
            y= .[[x]],
            split = .[['Farm']],
            type = 'violin',
            box=list(visible=T),
            meanline=list(visible=T)
          ) %>%
          layout(
            xaxis = list(
              title = "Population"
            ),
            yaxis = list(
              title = traits_title[[x]]
              #zeroline = F
            ),
            font=list(
              family = "sans serif",
              size = 20
            ),
            showlegend = FALSE
          )
        
        return(p)
      })
    }
  })
  
  size_plots<-reactive({
    req(input$traits)
    waiter<-Waiter$new(id='size_plot',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    selected_trait<-traits[traits_name%in%input$traits]
    if(input$population!='ALL'){
      req(input$fix)
      map2(selected_trait, input$fix,function(x,y){
        count_dat<-phe_farm() %>% 
          filter(!is.na(.data[[x]])) %>%
          group_by(.data[[y]]) %>% 
          summarise(Size=n())
        count_dat[[y]]<-factor(count_dat[[y]])
        p<-count_dat %>% 
          plot_ly(x=.[[y]],y=~Size,
                  color=.[[y]],
                  colors="#d62728",
                  opacity = 0.7) %>%
          add_bars() %>%
          layout(
            xaxis = list(
              title = y
            ),
            yaxis = list(
              title = paste0("Sample size")
              #zeroline = F
            ),
            font=list(
              family = "sans serif",
              size = 20
            ),
            showlegend = FALSE
          )
        return(p)
      })
    }else{
      map(selected_trait, function(x){
        count_dat<-phe_farm() %>% 
          filter(!is.na(.data[[x]])) %>%
          group_by(Farm) %>% 
          summarise(Size=n())
        count_dat[["Farm"]]<-factor(count_dat[["Farm"]],levels=farms)
        p<-count_dat %>% 
          plot_ly() %>%
          add_bars(x=~Farm,y=~Size,
                   color=~Farm,
                   colors=c("#d62728"),
                   opacity = 0.7) %>%
          layout(
            xaxis = list(
              title = "Line"
            ),
            yaxis = list(
              title = paste0("Sample size of ", x, " trait")
            ),
            font=list(
              family = "sans serif",
              size = 20
            ),
            showlegend = FALSE
          )
        
        return(p)
      })
    }
  })
  
  output$trait_plot<-renderPlotly(dens_plots()[[1]])
  output$size_plot<-renderPlotly(size_plots()[[1]])
}


shinyApp(ui = ui,server = server)

