# PHoloDB



## Introduction

Welcome to Pig hologenomic Data Base ([PHoloDB](http://alphaindex.zju.edu.cn/PHoloDB/)) v1.0!



[Hologenome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0457-9) means the sum of host genome and associated microbiome. It is a revolutionary concept proposed over ten years ago, and has provided a novel viewpoint to understand the formation, development, evolution, and characteristics of distinct biological entity. Pig is one of the most important meat suppliers and model animals for human beings. Therefore, it will be of great importance to reveal the genetic mechanisms of pigs' phenotypes from hologenomic perspective.



In order to help global scientists to make use of more integrated data to study the chracteristics of pigs using hologenomic methods and validate novel methods that can be used in pig breeding or genetic mechanism dissection, we develop PHoloDB to curate data, inluding but not limited to pig and its carried microorganism's genome, transcriptome, epigenome, phenotype, microbiome, metabolome, etc.



In the present version (1.0) of PHoloDB, we include the phenotype, genotype and microbiome data from three different pig lines. In the next version, we will further enlarge the sample size and also add the uploading function to make PHoloDB become a database that help global scientists to curate the pig hologenomics associated data.



## Function

Based on the curated data of version 1.0, PHoloDB at present include the following functions:



- Visualization of different omics data
- Data download and the explanation of the data structure can be found below



## Framework and source code

PHoloDB is nearly 100% developed by [Shiny](https://cran.r-project.org/web/packages/shiny/index.html) package in R and its friend packages. The source code of PHoloDB can be found within this repository.



## Installation

You can download the repository using the following command:

```shell
git clone git@github.com:ZheZhang-ZZ/PHoloDB.git
```



To run the app from your own computer, you need to install [R](https://cran.r-project.org/) first, and then run the following code within the R console to install the dependencies:

```R
packages_need_to_install<-c('shiny','shinydashboard','plotly','tidyverse','ggpubr','shinyjs','shinyWidgets','data.table','waiter','ggsci','leaflet','leaflet.extras','shinythemes','shinydashboardPlus')
packages_need_to_install <- packages_need_to_install[!(packages_need_to_install %in% installed.packages()[,"Package"])]
install.packages(packages_need_to_install)
```



If you have any trouble during installation, please propose an issue or contact me directly.



## Run app

We highly recommend you run the app within the [Rstudio](https://www.rstudio.com/) IDE, in which you can open the `app.R` script and just click the **Run App** button on the top right of the Rstudio window to get the code run in your browser or in a new window.



## Hologenomic data structure

The data can be downloaded within the [data base](http://alphaindex.zju.edu.cn/PHoloDB/) directly. In this repository, we include some example data which is in the same format with those in the PHoloDB.



### Phenotypic data

The phenotypic data files contain the following columns:

+ ID: The animal ID that is used through all dataset

+ Parity: The birth parity of the individuals

+ Herd: The herd in which the individuals are reared within the pig farm, please note that the herd is recoded using natural numbers, so the meaning of these numbers only maintain within each file

+ YearSeason: A number with three digit, in which the first two digits are the year and the last digit is the season (1- sprint; 2- summer; 3- autumn; 4- winter)

+ Sex: The sex of the individuals

+ Dam: The dam of each individual

+ Birth_weight: The birth weight of pigs, which were recorded within 24 hours after the piglets were born

+ Day_to_115kg: The number of days to reach 115 kg, which was calucated based on the body weight at the off-test stage using the following formula:

  
  $$
  \text {Adjusted days to 115 kg} = \text {actual age} + \left[ 115 - \text{actual weight} \times \frac{\text{actual age}-a}{\text{actual weight}} \right]
  $$
  

  where $a=50$ for boars and $40$ for sows

+ Backfat_thickness: The backfat thickness (mm) of pigs were measured using ultrasound scanning between the 3rd and 4th last ribs, which were corrected to 115 kg using the following formula:

  
  $$
  \text{Adjusted backfat} = \text{actual backfat} + \left[\left(115-\text{actual weight}\right)\times\frac{\text{actual backfat}}{\text{actual weight}-a}\right]
  $$
  

  where $a=-9.07$ for boars and $2.68$ for sows

+ Loin_muscle_area: The loin muscle area ($\text{cm}^2$) of pigs were measured  using ultrasound scanning between the 3rd and 4th last ribs, which were corrected to 115 kg using the following formula:

  
  $$
  \begin{aligned}
  \text{Adjusted loin muscle area} &= \text{actual loin muscle area} +\\
  & \left[\left(115-\text{actual weight}\right)\times\frac{\text{actual loin muscle area}}{\text{actual weight}+ 70.31}\right]
  \end{aligned}
  $$
  

+ Lean_percentage: The  lean percentage (%) were predicted based on backfat thickness and loin muscle area ajusted to 115 kg using the following formula:

  
  $$
  \begin{aligned}
  \text {lean percentage} &= \left[80.96-\left(0.6472 \times \text{Adjusted backfat}\right) + \\
  \left( 0.7274 \times \text{Adjusted loin muscle area}\right)\right] \times 0.54
  \end{aligned}
  $$



### Pedigree data

The pedigree file contains three columns:

+ ped_ind: The individual ID
+ sire_ind: The sire of the individual
+ dam_ind: The dam of the individual



### Genotypic data

The raw fastq files of sequence data are available at  SRA database under ID SRP323506 within the BioProject [PRJNA736175](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA736175/). We include tow sets of genotypic data here. One is the SNP array dataset in which all chromosome loci are contained in Plink [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed) format, and the other one is the imputed SNP dataset which contains15142178 autosomal SNPs in different chromosome-wise files in Plink [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed) format.



### Microbiomic data

The raw fastq files of shotgun metagenomic data are available at SRA database under ID SRP323499 within the BioProject [PRJNA736198](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA736198). The absolute taxonomic counts calculated by Kraken2 and relative abundances (%) profiled by MetaPhlAn3 at different levels (phylum, class, order, family, genus and species) are available in PHoloDB. 

The count data are in files with name `*count*.txt` while abundance data are in files with names `*abundance*.txt`. The two files' data structure are self-explanatory, in which the columns are taxon names at different levels and individual ID.



## Citation

We hope that the citation paper will be coming soon...

