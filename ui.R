# MOPower v1.0.2
# User Interface
# You can run the application by clicking 'Run App' above.

library(CancerSubtypes)
library(extRemes)
library(ggplot2)
library(glmpath)
library(lme4)
library(mediation)
library(NMF)
library(omicade4)
library(plotly)
library(reshape)
library(reshape2)
library(shiny)
library(shinyBS)
library(shinythemes)
library(shinydashboard)
library(SNFtool)
library(MOFA)
library(survival)
library(reticulate)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(iClusterPlus)
library(iCluster)
library(BiocInstaller)
library(BiocVersion)


## selection of omics tickboxes
# Define UI for application
shinyUI(fluidPage(
  # Application title
  titlePanel(tags$head(tags$link(rel = "icon", type = "image/png", href = "MOPowerlogo4.png"),
                       tags$title("MOPower"))
  ),
  navbarPage(theme = shinythemes::shinytheme("spacelab"),  title=div(img(src="MOPowerlogo4.png", height = 50), "MOPower"),

  tabPanel("PC", h2("Multi-Omics Power Calculator & Data simulator"), 
  fluidRow(
     column(12,wellPanel(
               helpText("This power calculator is used for the design of multi-omic studies. 
                        Simulation of genomic, transcriptomic and epigenomic data with case-control 
                        or time-to-event outcomes. Power calculation is based on the power to detect 
                        a locus association incorporating multi-omic festures with the outcome of interest. 
                        Single omics power and sample size calculation can also be undertaken using MOPower."),
             
               h4("Select and fill in parameter values"),
             # Selecting omics
               checkboxGroupInput(inputId = "omics", label ="Omics selection", choices = list("Genome", "RNA-seq", "Epigenome"
                                                                                              #,"Proteome","Metabolome"
                                                                                              ),inline=TRUE),
               bsTooltip("omics", "Genomic data - SNPs, Rare variants & CNVs. \\ RNA-Seq data - Gene expression (Read counts). \\ Epigenome data - DNA methylation (array). \\ Proteomic data - Protein counts.","bottom", options = list(container = "body")),
            # Selecting outcome
            checkboxGroupInput(inputId = "outselect", label ="Outcome/Phenotype selection", choices = list("Case-Control (Binary)"= 1, "Time-to-event (Survival)"= 3),inline=TRUE),
            bsTooltip("outselect", "Choose a single outcome","bottom", options = list(container = "body"))
            
            ))
     ),
  # Conditional reactive panels
  fluidRow(
     column(2,style='padding:0px;',wellPanel(
            h4("Study Design Parameters"),
            conditionalPanel(
              condition = "(input.omics[0] != 'RNA-seq' || input.omics[1] == 'RNA-seq') && (input.select != 19 || input.select != 20)",
              numericInput("replicates", "Simulations",value = 1, min = 1, max = 10000)
            ),
            numericInput("obs", "Sample size", min = 1, max = 1000, value = 40),
            conditionalPanel(
            condition = "input.outselect == 1",
            numericInput("ccratio", "Case-Control Ratio", value = 0.5, min = 0, max = 1)
            ),
            conditionalPanel(
            condition = "input.outselect == 3",
            selectInput("timedist", label = "Survival Distribution", choices = list("Weibull" = 1, "Log-Normal" = 2, "Gamma" = 3), selected = 1),
            numericInput("Stime", "Survival Time (Baseline)", value = 20, min = 0, max = 1000),
            selectInput("censtype", label = "Type of Censoring", choices = list("Right" = 1
                                                                                #, "Left" = 2, "Interval" = 3
                                                                                ), selected = 1),
            numericInput("Ctime", "Percentage of censored Obs.", value = 20, min = 0, max = 1000),
            bsTooltip("Ctime", "This is highly dependent on the baseline hazard and end of study time.","right", options = list(container = "body")),
            numericInput("eos", "End of study (Days)", value = 100, min = 0, max = 1000)
            ),
            numericInput("tratio", "Treatment Allocation Ratio", value = 0, min = 0, max = 1),
            numericInput("treat", "Treatment Effect Size", value = 0, min = 0, max = 40),
          
       conditionalPanel(
         condition = "input.outselect == 2",
            numericInput("fusamples", "Follow-up samples taken (Every __ Days)", value = 50, min = 0, max = 100)
       )
     )),
     column(2,wellPanel(h4("Gene-level data simulation Parameters"),
      conditionalPanel(
      condition = "input.omics[0] == 'Genome' || input.omics[1] == 'Genome' || input.omics[2] == 'Genome'",
      sliderInput("snps", "Variants", min = 1, max = 1000, value = 4),
      sliderInput("snpeffect", "SNP effect range", min = 0.00001, max = 10.0000, value = c(0.1,0.2)),
      #numericInput("snpeffect", "Average SNP Effect Size", value = 0.4, min = 0, max = 40),
      sliderInput("mafint", "MAF range", min = 0.000001, max = 0.5, value = c(0.01,0.2))
      ),
   #  h4("Gene-level data simulation Parameters"),
    conditionalPanel(
      condition = "input.omics[0] == 'RNA-seq' || input.omics[1] == 'RNA-seq' || input.omics[2] == 'RNA-seq'",
      numericInput("genes", "Genes", min = 1, max = 100000, value = 5),
      numericInput("reads", "Average baseline reads per gene", value = 100, min = 0, max = 500),
      sliderInput("fc", "Fold Change (log2)", min = 0.1, max = 5.0, value = c(0.5,1.5)),
      sliderInput("disp", "Dispersion", min = 0.01, max = 5.00, value = c(0.5,2.5)),
     # numericInput("fc", "Fold Change (log2)", value = 2, min = -100, max = 100),
     # numericInput("disp", "Dispersion", value = 0.5, min = 0, max = 100),
      bsTooltip("disp", "Usual values for the square root dispersion (D) are between 0.1 and 1. The dispersion equation used in the NB distribution is 1/(D) ","right", options = list(container = "body"))
    ),
  #  h4("CpG site data simulation Parameters"),
    conditionalPanel(
      condition = "input.omics[0] == 'Epigenome' || input.omics[1] == 'Epigenome' || input.omics[2] == 'Epigenome'",
      sliderInput("cpg", "Number of CpG sites", min = 1, max = 10000, value = 5),  
      numericInput("casemeth", "Methylation beta distribution values for cases", min = 1, max = 10, value = 6),
      numericInput("casemeth2","", min = 1, max = 10, value = 2),
      numericInput("contmeth", "Methylation beta distribution values for controls", min = 1, max = 10, value = 2),
      numericInput("contmeth2","", min = 1, max = 10, value = 6),
      textOutput("meanmeth"),
      checkboxGroupInput("checkmeth","", choices = list("Unmethylated sites" = 1, "Hemi-Methylated sites" = 2))
    ),
    conditionalPanel(
      condition = "input.checkmeth[0]",
      sliderInput("unmeth", "Number of Un-methylated sites", min = 1, max = 10000, value = 1)
    ),
    conditionalPanel(
      condition = "input.checkmeth[1]",
      sliderInput("hemimeth", "Number of Hemi-methylated sites", min = 1, max = 10000, value = 1)
    )
    )),
    
    # Methods selection
    column(3,wellPanel(
      selectInput("select", h4("Statistical integration/analysis models"), choices = c()),
      bsTooltip("select", "Combination of multi-omics data integration methods & regression-based models for analysis","right", options = list(container = "body"))
      
      )),
    # Variable selection for analysis models
    column(2, checkboxGroupInput("checkGroup",  h5("Parameters"), 
                       choices = c()),
    
    conditionalPanel(
             condition = "input.select == 5 || input.select == 17",
    selectInput("mediator","Select mediator variable", choices = c())
    ),
    conditionalPanel(
      condition = "input.select == 5 || input.select == 17",
      selectInput("causalcov","Select dependent causal variable", choices = c())
    ),
    conditionalPanel(
      condition = "input.select == 2 || input.select == 3",
      numericInput("mofathresh","Variance explained threshold for factors", value = 0.02, min = 0, max = 1),
      bsTooltip("mofathresh", "e.g. 0.02 = 2% of variance explained. All factors contributing less than 0.02 are removed from the analysis.","right", options = list(container = "body"))
    ),
    numericInput("error", "Significance threshold (Type I error)", value = 0.01, min = 0, max = 1)
    ),
    column(2,
    downloadButton("dloadData", "Download Sample Data"),
    column(2,br(),
#downloadButton("report", "Generate report"), br(), br(), 
actionButton("submit","Power Calculation", style = "color: white; 
                       background-color: green; height: 60px; width: 150px; text-align:center; text-indent: -2px; border-radius: 6px"))
    ),
       mainPanel(tags$style(type="text/css",
                            ".shiny-output-error { visibility: hidden; }",
                            ".shiny-output-error:before { visibility: hidden; }"
       ),
         tabsetPanel(
           tabPanel("Data Replicate", tableOutput('table')), 
         # tabPanel("Replicate Analysis Summary", verbatimTextOutput("summary")),
           tabPanel("Power Summary", verbatimTextOutput("powersummary")), 
           tabPanel("Plots", fluidRow(br(),
                    column(3,textInput("label", "Variable label", "")),
                    checkboxInput("addline", "Add line", FALSE),
                    column(12,plotlyOutput("plot")))
                    )
         )
       )
       )
       ),
  # False discovery rate page
      tabPanel("FDR & Bonferroni", fluidRow(
        column(8,wellPanel(h5("In sequencing experiments we are fitting one model for each gene/exon/sequence of interest, and therefore performing thousands
               of tests. We need to correct for multiple testing using the Bonferroni correction and calculate the FDR.")))),
               fluidRow(
                 column(12,mainPanel(
                   h3("Bonferroni correction"),
                   numericInput("genes2", "Total number of genes, SNPs, CpG sites or features", min = 1, max = 1000000, value = 100),
                   numericInput("error2", "Significance threshold", value = 0.01, min = 0, max = 1),
                   h6("Bonferroni correction for number of hypotheses tested"),
                   verbatimTextOutput("BC"),
                   h3("Distribution of calculated p-values"),
                   plotlyOutput("distPlot"),
                 conditionalPanel(
                 condition = "input.select == 19 || input.select == 20",
                 h3("Estimated dispersion"),
                 plotlyOutput("disphist")
                   )
                 )), column(12,mainPanel(
                   h3("FDR (False discovery rate)"),
                   verbatimTextOutput("false"),
                   plotlyOutput("fdrplot"))
                 )
                 )
                 
                 ),
    #  tabPanel("Cost Analysis", "Cost Analysis based on simulated study design."),
      #tabPanel("Instructions", fluidRow(
        #column(10,wellPanel(h3("Example"),h4("We are designing a time-to-event study where we follow-up 20 patients for 100 days after recruitment. 
                   #                         We want to test the power to detect an association with a specific gene that countains 15 SNPs and 
                  #                          10 CpG sites with a 0.071 difference in mean methylation between patients who experience the event and those that are censored. We have gene expression data collected for each patient. Using MOPower we 
                 #                           can calculate power under different scenarios, changing the MAF interval for the SNPs, 
                #                            adding unmethylated CpG sites and testing multiple integration methods, however for the purpose of this example we will just show 
               #                             one power calculation using the Cox Path model with L1 Penalty.")),
              # mainPanel(img(src="1.png",height = 800),img(src="4.png", height = 400),img(src="5.png", height = 800))))),
      tabPanel("About", 
               fluidRow(
                 column(8,wellPanel(h3("MOPower v1.02, GNU General Public Licence"),"MOPower is a multi-omics power calculator, implemented with the latest integration and analysis models. Data is simulated using statistical distributions and user input parameters. Power is calculated based on the number of replicates
               with a significant association for total effects of omics datasets with the outcome of interest (case-control or survival). See vignette for more details: https://github.com/HSyed91/MOPower. MOPower was developed at GOSgene, Institute of Child Health, University College London. This research is funded by the BRC."
                 ), img(src="NIHR.jpg",height = 50),"----",img(src="GOSgene.png", height = 50)
                 ))),
      tabPanel("News", fluidRow(
        column(8,wellPanel(h3("Updates, new features & bug fixes")),
               h5("06/08/2019 v1.02"),"1. Parameters such as effect sizes have been changed from single value for all simulated features to range values for a more realistic setting.",br(),br(),
               "2. Multiprocessing using 5 CPU cores has been added for concurrent power calculation.", br(),br(),
               h5("10/04/2019 v1.02"),"1. Single omics power calculation of RNA-Seq data can now be calculated using edgeR's Exact test or Negative Binomial GLM. The output is now the power to detect differentially expressed genes between groups given a sample of genes (Q: How many genes will be differentially expressed?) Number of simulations input will disappear as this is not required.",br(),br(),
               "2. Plot of estimated dispersion parameter is displayed for single omics RNA-Seq power calculation.", br(),br(), 
               "3. An 'Add line' checkbox to add a new power calculation line to the power plot has been added.",br(),br(),
               "4. Data replicate displayed in UI will now only show one feature from each omics. Download sample data will have the entire simulated dataset. This was too minimize the space taken by the display data.",br(),br(),
               "5. Longitudinal outcome removed from simulation framework because there are no multi-omics analyses currently available to analyse the outcome.",br(),br(),
               "6. Analysis options shown in menu are dependent on the outcome choice.",br(),br(),
               "7. Dispersion parameter within the Negative-Binomial model used to simulate gene expression data is now defined using the edgeR definition of 1/(D^2)"))),
      tabPanel("FAQ", fluidRow(
        column(8,wellPanel(h3("Frequently Asked Questions"), "Please submit your questions via email to h.syed@ucl.ac.uk. Questions and answers will be posted on this page once the question has been answered or the issue has been resolved.")))),
      tabPanel("References", fluidRow(
        column(8,wellPanel(h3("References"),
"1. Ching, T., Huang, S., & Garmire, L. X. (2014). Power analysis and sample size estimation for RNA-Seq differential expression. RNA, 20(11), 1684-1696. doi:10.1261/rna.046011.114.",br(),

"2. Li, C. I., & Shyr, Y. (2016). Sample size calculation based on generalized linear models for differential expression analysis in RNA-seq data. Stat Appl Genet Mol Biol, 15(6), 491-505. doi:10.1515/sagmb-2016-0008",br(),

"3. Miller, F., Zohar, S., Stallard, N., Madan, J., Posch, M., Hee, S. W., . . . Day, S. (2018). Approaches to sample size calculation for clinical trials in rare diseases. Pharm Stat, 17(3), 214-230. doi:10.1002/pst.1848 ",br(),

"4. Plan, E. L. (2014). Modeling and simulation of count data. CPT Pharmacometrics Syst Pharmacol, 3, e129. doi:10.1038/psp.2014.27",br(),

"5. Poplawski, A., & Binder, H. (2017). Feasibility of sample size calculation for RNA-seq studies. Brief Bioinform. doi:10.1093/bib/bbw144",br(),

"6. Syed, H., Jorgensen, A. L., & Morris, A. P. (2016). SurvivalGWAS_Power: a user friendly tool for power calculations in pharmacogenetic studies with 'time to event' outcomes. BMC Bioinformatics, 17(1), 523. doi:10.1186/s12859-016-1407-9",br(),

"7. Tsai, P. C., & Bell, J. T. (2015). Power and sample size estimation for epigenome-wide association scans to detect differential DNA methylation. Int J Epidemiol, 44(4), 1429-1441. doi:10.1093/ije/dyv041",br(),

"8. Wu, H., Wang, C., & Wu, Z. (2015). PROPER: comprehensive power evaluation for differential expression using RNA-seq. Bioinformatics, 31(2), 233-241. doi:10.1093/bioinformatics/btu640",br(),

"9. Yu, L., Fernandez, S., & Brock, G. (2017). Power analysis for RNA-Seq differential expression studies. BMC Bioinformatics, 18(1), 234. doi:10.1186/s12859-017-1648-2",br(),

"10. Argelaguet, R., Velten, B., Arnol, D., Dietrich, S., Zenz, T., Marioni, J. C., . . . Stegle, O. (2018). Multi-Omics Factor Analysis-a framework for unsupervised integration of multi-omics data sets. Mol Syst Biol, 14(6), e8124.",br(),

"11. Huang, S., Chaudhary, K., & Garmire, L. X. (2017). More Is Better: Recent Progress in Multi-Omics Data Integration Methods. Front Genet, 8, 84. doi:10.3389/fgene.2017.00084",br(),

"12. Rohart, F., Gautier, B., Singh, A., & Lê Cao, K. A. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. PLoS Comput Biol, 13(11), e1005752. doi:10.1371/journal.pcbi.1005752 ",br(),

"13. Sun, Y. V., & Hu, Y. J. (2016). Integrative Analysis of Multi-omics Data for Discovery and Functional Studies of Complex Human Diseases. Adv Genet, 93, 147-190. doi:10.1016/bs.adgen.2015.11.004",br(),

"14. Wang, B., Mezlini, A. M., Demir, F., Fiume, M., Tu, Z., Brudno, M., . . . Goldenberg, A. (2014). Similarity network fusion for aggregating data types on a genomic scale. Nat Methods, 11(3), 333-337. doi:10.1038/nmeth.2810",br(),

"15. https://www.bioconductor.org/packages/devel/bioc/vignettes/CancerSubtypes/inst/doc/CancerSubtypes-vignette.html"))))
     )
   )
  )



