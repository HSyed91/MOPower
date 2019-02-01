# MOPower v1.0.1
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.

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
library(SNFtool)
library(survival)
library(reticulate)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(iClusterPlus)
library(iCluster)
library(BiocInstaller)
library(BiocVersion)
library(MOFAtools)

# v This section is for shinyapps.io v
#virtualenv_create("r-reticulate")
#use_virtualenv("r-reticulate")
#use_python("/usr/bin/python3.5")
#py_install("mofapy")

# Define server logic
shinyServer(function(input, output, session) {
  # Bonferroni correction value
  output$BC <- renderPrint(
    {
      input$error2/input$genes2
    }
  )
  datainputs <- reactive({
  # Study design parameters
  # Reactive datatable
  n <- input$obs
  EndofStudy <- input$eos
  Treatment <- rbinom(n,1,input$tratio)
  # Gene-level parameters
  numsnp <- input$snps
  numgene <- input$genes
  Patients <- 1:n
  sum1 <- -0
  dt <- data.frame(Patients,Treatment,EndofStudy)
  if("Genome" %in% input$omics){
  # SNPs
  snpdt <- as.data.frame(matrix(NA, ncol = numsnp, nrow = n))
  beta <- c(1:numsnp)
  agg<- as.data.frame(matrix(NA, ncol = 4, nrow = n))
  linpred <- matrix(NA, nrow = n,ncol=numsnp)
  # SNP data
  for(i in 1:numsnp){
    maf <- runif(1,input$mafint[1], input$mafint[2])
    snpdt[i] <- cbind(rbinom(n,2,maf))
    colnames(snpdt)[i] <- paste("SNP", i, sep = "")
    beta[i] <- sqrt(input$snpeffect/2/maf/(1-maf))
  }
  agg <- t(t(snpdt)*beta)
  sum1 <- apply(agg, 1, sum)
  dt <- cbind(dt,snpdt)
  }
  # Outcome (case-control)
  if(input$outselect == 1){
    if("Genome" %in% input$omics){
  b0 <- log(input$ccratio/(1-input$ccratio)) -  sum1 - (Treatment*input$treat)
  linpred <- b0 + sum1 + (Treatment*input$treat)
  pi <- exp(linpred) / (1 + exp(linpred))
  Group <- rbinom(n, size=1, prob=pi)
    }
    else{
      Group <- rbinom(n, size=1, prob=input$ccratio)
    }
  dt$Group <- Group
  }
  # Outcome (Survival)
  else if(3 %in% input$outselect){
  # baseline hazard
  lambda <- input$Stime*exp(-sum1 - (Treatment*input$treat))
  survtime <- rweibull(input$obs,1, scale=lambda)
  censtime <- rbinom(input$obs,1,1-(input$Ctime*0.01))
  # Right censoring as percentage
  Group <- ifelse(survtime >= input$eos,0,censtime)
  Observedtime <- ifelse(survtime>=input$eos,input$eos,survtime)
  # v Censoring is random v
  #Observedtime <- ifelse(survtime>input$eos & censtime>input$eos,input$eos,pmin(survtime,censtime))
  #Group <- ifelse(survtime < censtime & survtime < input$eos,1,0)
  dt$Censoring <- Group
  dt$Observedtime <- Observedtime
  }
  # Outcome (longitudinal)
  else if(2 %in% input$outselect){
    timevec <- matrix(rep(seq(0, input$eos, input$fusamples), input$obs), ncol=(input$eos/input$fusamples)+1, nrow=input$obs, byrow=T)
    colnames(timevec) <- colnames(timevec, do.NULL = FALSE, prefix = "Time")
    dt <- cbind(dt,timevec)
  }
  # Gene expression
  if("RNA-seq" %in% input$omics){
  mexp <-  exp(Group*input$fc+log(input$reads)) # N can be input or estimated from uniform distribution. N = library size
  disp <- input$disp #  or variance = mexp*(1+mexp*input$disp)or variance=mexp+(mexp^2/input$disp)
  GE <- as.data.frame(matrix(NA, ncol = numgene, nrow = n))
  for(j in 1:numgene){
    GE[j] <- cbind(rnbinom(n, size=disp, mu=mexp))
    colnames(GE)[j] <- paste("GE", j, sep = "")
  }
  # v Time vector for gene expression over time v
  #timevec <- matrix(rep(seq(0, input$eos, input$fusamples), input$obs), ncol=(input$eos/input$fusamples)+1, nrow=input$obs, byrow=T)
  #colnames(timevec) <- colnames(timevec, do.NULL = FALSE, prefix = "Time")
  
  # DT without Methylation
  dt <- cbind(dt,GE)
  }
  
  # DNA methylation - array CpG sites
  if("Epigenome" %in% input$omics){
  sites <- input$cpg
  # r and v are unmeth and hemimeth are column numbers
  r <- 0
  v <- 0
  unsites <- as.data.frame(matrix(NA, ncol = input$unmeth, nrow = n))
  hemisites <-as.data.frame(matrix(NA, ncol = input$hemimeth, nrow = n))
  methsites <- as.data.frame(matrix(NA, ncol = sites, nrow = n))
  for(k in 1:sites){
    meth = c(1:n)
    for(m in 1:length(Group)){ 
      if(Group[m]==1){
        meth[m] <- rbeta(1,input$casemeth,input$casemeth2)
      }
      else{
        meth[m] <- rbeta(1,input$contmeth,input$contmeth2)
      }
    }
    methsites[,k] = cbind(meth)
    colnames(methsites)[k] <- paste("CpGmeth", k, sep = "")
  }
 # unmeth and hemi meth sites
  if(1 %in% input$checkmeth){
  for(v in 1:input$unmeth){
  unsites[,v] <- rbeta(n,1,6)
  colnames(unsites)[v] <- paste("CpGmeth", v+k+r, sep = "")
  }
    methsites <- cbind(methsites,unsites)
  }
    if(2 %in% input$checkmeth){
  for(r in 1:input$hemimeth){
    hemisites[,r] <- rbeta(n,6,6)
    colnames(hemisites)[r] <- paste("CpGmeth", r+v+k, sep = "")
  }
      methsites <- cbind(methsites,hemisites)
    }

  dt <- cbind(dt,methsites)
  }

  # Print sample data
   dt
  })
  output$table <- renderTable({
    
    head(datainputs(), n = 20)
    })   

 #power calc button 
 #reactive({
 #observeEvent(input$submit,{session$sendCustomMessage(type = 'message', message = 'Power is being calculated. Please wait....')
 # })
 #}) 
 
  # Integration and analysis models
  methods <- eventReactive(input$submit, {
  pcaldf <- as.data.frame(matrix(NA, ncol = 2, nrow = 10))
  pdf2 <- as.data.frame(matrix(NA, ncol = 2, nrow = input$replicates*10))
  pdfadd <- as.data.frame(matrix(NA, ncol = 2, nrow = 10))
  
  withProgress(message = 'Calculating Power', value = 0, {
  for (z in 1:10){ # power at 10 intervals on plot
  value <- 0
  value2 <- 0
  pdf <- as.data.frame(matrix(NA, ncol = 1, nrow = input$replicates))
  pd <- as.data.frame(matrix(NA, ncol = 1, nrow = input$replicates))
  
  #Simulate all data replicates  
  for (h in 1:input$replicates){
 
  n <- input$obs*z
  sites <- 0
  unsites <- 0
  hemisites <- 0
  EndofStudy <- input$eos
  methset <- NULL
  geneset <- NULL
  snpset <- NULL
  sum1 <- 0
  Treatment <- rbinom(n,1,input$tratio)
  # Gene-level parameters
  numsnp <- input$snps
  numgene <- input$genes
      Patients <- 1:n
      
      if("Genome" %in% input$omics){
      #SNPs
      snpdt <- as.data.frame(matrix(NA, ncol = numsnp, nrow = n))
      beta <- c(1:numsnp)
      agg<- as.data.frame(matrix(NA, ncol = 4, nrow = n))
      linpred <- matrix(NA, nrow = n,ncol=numsnp)
      for(i in 1:numsnp){
        maf <- runif(1,input$mafint[1], input$mafint[2])
        snpdt[i] <- cbind(rbinom(n,2,maf))
        colnames(snpdt)[i] <- paste("SNP", i, sep = "")
        beta[i] <- sqrt(input$snpeffect/2/maf/(1-maf))
      }
      agg <- t(t(snpdt)*beta)
      sum1 <- apply(agg, 1, sum)
      dt <- data.frame(Patients,Treatment,EndofStudy,snpdt)
      allsnp <- dt[ , grepl( "SNP" , names( dt ) ) ]
      snpset <- apply(allsnp, 1, sum)
      }
        
      # Outcome (case-control)
      if(input$outselect == 1){
        if("Genome" %in% input$omics){
        b0 <- log(input$ccratio/(1-input$ccratio)) -  sum1 - (Treatment*input$treat)
        linpred <- b0 + sum1 + (Treatment*input$treat)
        pi <- exp(linpred) / (1 + exp(linpred))
        Group <- rbinom(n, size=1, prob=pi)
      }
      else{
        Group <- rbinom(n, size=1, prob=input$ccratio)
      }
        dt$Group <- Group
      }
      
      # Outcome (survival)
      else if(3 %in% input$outselect){
        # baseline hazard
        lambda <- input$Stime*exp(-sum1 - (Treatment*input$treat))
        survtime <- rweibull(input$obs*z,1, scale=lambda)
        censtime <- rbinom(n,1,1-(input$Ctime*0.01))
        # right censoring
        Group <- ifelse(survtime >= input$eos,0,censtime)
        Observedtime <- ifelse(survtime>=input$eos,input$eos,survtime)
        dt$Group <- Group
        dt$Observedtime <- Observedtime
      }
      
      # Gene expression
      if("RNA-seq" %in% input$omics){
      mexp <-  exp(Group*input$fc+log(input$reads)) 
      disp <- input$disp 
      GE <- as.data.frame(matrix(NA, ncol = numgene, nrow = n))
      for(j in 1:numgene){
        GE[j] <- cbind(rnbinom(n, size=disp, mu=mexp))
        colnames(GE)[j] <- paste("Genecount", j, sep = "")
      }
      dt <- cbind(dt,GE)
      allGE <- dt[ , grepl( "Genecount" , names( dt ) ) ]
      geneset <- apply(allGE, 1, sum)
      }
     # timevec <- matrix(rep(seq(0, input$eos, input$fusamples), n), ncol=(input$eos/input$fusamples)+1, nrow=n, byrow=T)
     # colnames(timevec) <- colnames(timevec, do.NULL = FALSE, prefix = "Time")
      
      # DNA methylation - array CpG sites
      if("Epigenome" %in% input$omics){
        sites <- input$cpg
        v <- 0
        r <- 0
        unsites <- as.data.frame(matrix(NA, ncol = input$unmeth, nrow = n))
        hemisites <-as.data.frame(matrix(NA, ncol = input$hemimeth, nrow = n))
        methsites <- as.data.frame(matrix(NA, ncol = sites, nrow = n))
        for(k in 1:sites){
          meth = c(1:n)
          for(m in 1:length(Group)){ 
            if(Group[m]==1){
              meth[m] <- rbeta(1,input$casemeth,input$casemeth2)
            }
            else{
              meth[m] <- rbeta(1,input$contmeth,input$contmeth2)
            }
          }
          methsites[,k] = cbind(meth)
          colnames(methsites)[k] <- paste("CpGmeth", k, sep = "")
        }
        # unmeth and hemi meth sites
        if(1 %in% input$checkmeth){

          for(v in 1:input$unmeth){
            unsites[,v] <- rbeta(n,1,6)
            colnames(unsites)[v] <- paste("CpGmeth", v+k+r, sep = "")
          }
          methsites <- cbind(methsites,unsites)
        }
          if(2 %in% input$checkmeth){
          for(r in 1:input$hemimeth){
            hemisites[,r] <- rbeta(n,6,6)
            colnames(hemisites)[r] <- paste("CpGmeth", r+v+k, sep = "")
          }
            methsites <- cbind(methsites,hemisites)
          }

        
        dt <- cbind(dt,methsites)
        allmeth <- dt[ , grepl( "CpGmeth" , names( dt ) ) ]
        methset <- apply(allmeth, 1, sum)
      }
      # Joining data
    
    # Progress bar
    incProgress(1/z, detail = paste("Please wait...."))
  
  # Methods
  # Joint regression model
  if(input$select==1){
    dta <- data.frame(Group)
    if (is.null(geneset) == FALSE){
      dta <- cbind(dta,geneset)
    }
    if (is.null(snpset)== FALSE){
      dta <- cbind(dta,snpset)
    }
    if (is.null(methset)== FALSE){
      dta <- cbind(dta,methset)
    }
    if (input$tratio > 0 && 4 %in% input$checkGroup){
      dta <- cbind(dta,Treatment)
    }
 
  # Logistic regression
    model <- glm(Group~.,data=dta,family=binomial(link='logit'))
    sumup <- summary(model)
    model2 <- glm(Group~1,data=dta,family=binomial(link='logit'))
    anova <- anova(model,model2,test="Chisq")  
    pval1 <- anova[2,5]
    pdf[h,] <- pval1
 #  pdf2[h,2] <- pval1

    if(pdf[h,]<=input$error){
      value <- value + 1
    }
  }
    
    if(input$select==16){
      dta <- data.frame(Group)
      if (is.null(geneset) == FALSE){
        dta <- cbind(dta,geneset)
      }
      if (is.null(snpset)== FALSE){
        dta <- cbind(dta,snpset)
      }
      if (is.null(methset)== FALSE){
        dta <- cbind(dta,methset)
      }
      if (input$tratio > 0 && 4 %in% input$checkGroup){
        dta <- cbind(dta,Treatment)
      }
      
      # Cox PH regression
      model <- coxph(Surv(Observedtime, Group) ~.,data=dta)
      sumup <- summary(model)
      model2 <- coxph(Surv(Observedtime, Group) ~1,data=dta)
      anova <- anova(model,model2,test="Chisq")  
      pval1 <- anova[2,5]
      pdf[h,] <- pval1
      if(pdf[h,]<=input$error){
        value <- value + 1
      }
    }
    
  # Cox Path L1 Penalty
    if(input$select==10){
      
      features <- data.frame(matrix(ncol=0,nrow=n))
      if (is.null(geneset) == FALSE){
        features <- cbind(geneset)
      }
      if (is.null(snpset)== FALSE){
        features <- cbind(features,snpset)
      }
      if (is.null(methset)== FALSE){
        features <- cbind(features,methset)
      }
      if (input$tratio > 0 && 4 %in% input$checkGroup){
        features <- cbind(features,Treatment)
      }
      dat <- list(x=features, time=Observedtime, status=Group)
      invisible(capture.output(fit.a <- coxpath(dat)))
      len <- length(fit.a$loglik)
      last <- tail(fit.a$df, n=1)
      b <- lr.test(fit.a$loglik[len], fit.a$loglik[1], alpha = input$error, df = last)
      pdf[h,] <- b$p.value
      if(z==10 && h==input$replicates){
       sumup <- summary(fit.a) 
      }
      if(pdf[h,]<=input$error){
      value <- value + 1
      } 
    }
    
    # mediation analysis
    if(input$select==5 || input$select == 17){
      
      if(2 %in% input$mediator){
        mediator <- geneset
      }
      else if(1 %in% input$mediator){
        mediator <- snpset
      }
      else {
        mediator <- methset
      }
      
      if(2 %in% input$causalcov){
        abc <- geneset
      }
      else if(1 %in% input$causalcov){
        abc <- snpset
      }
      else {
        abc <- methset
      }
      
      dta <- data.frame(Group,mediator,abc)
      
      if(1 %in% input$checkGroup && ((2 %in% input$mediator && 3 %in% input$causalcov) || (3 %in% input$mediator && 2 %in% input$causalcov))){
        b <- snpset
        dta <- cbind(dta,b)
        }
      else if(2 %in% input$checkGroup && ((1 %in% input$mediator && 3 %in% input$causalcov)|| (3 %in% input$mediator && 1 %in% input$causalcov))){
        b <- geneset
        dta <- cbind(dta,b)
      }
      else if(3 %in% input$checkGroup && ((2 %in% input$mediator && 1 %in% input$causalcov) || (1 %in% input$mediator && 2 %in% input$causalcov))){
        b <- methset
        dta <- cbind(dta,b)
      }

      if (4 %in% input$checkGroup){
      c <- Treatment
      dta <- cbind(dta,c)
      }
      
      med.fit <- lm(mediator ~.,data=dta[,!(names(dta) %in% c('Group'))])
      if(input$select==5){
      out.fit <- glm(Group~.*., data=dta,family = binomial("probit"))
      }
      else{
      out.fit <- survreg(Surv(Observedtime,Group)~.*., data = dta, dist='weibull')
      }
      try(
      med.out <- mediate(med.fit, out.fit, treat = "abc" , mediator = "mediator",robustSE = TRUE, sims = 100)
      )
      sumup <- summary(med.out)
      a <- summary(med.out)
      # total effect p-value
      pdf[h,] <- a$tau.p
      if(pdf[h,]<=input$error){
        value <- value + 1
      } 
      
    }
    
    #mofa
   if(input$select==2 || input$select==3){
  
     a <- list()
     
     if (is.null(geneset) == FALSE){
       if(h==1 && z==1){
         rownames(allGE) <- paste("Sample", 1:input$obs*z, sep = "")
       }
       a <- list(t(allGE))
       names(a)<-c("GE")
     }
     if (is.null(snpset)== FALSE){
       if(h==1 && z==1){
         rownames(allsnp) <- paste("Sample", 1:input$obs*z, sep = "")
       }
       a <- append(a,list(t(allsnp)))
       
       if (is.null(geneset) == FALSE){
       names(a)<-c("GE","SNP")
       }
       else{
         names(a)<-c("SNP")
       }
     }
     if (is.null(methset)== FALSE){
       if(h==1 && z==1){
         rownames(allmeth) <- paste("Sample", 1:input$obs*z, sep = "")
       }
       a <- append(a,list(t(allmeth)))
       
       if (is.null(geneset) == FALSE && is.null(snpset) == FALSE){
         names(a)<-c("GE", "SNP", "Meth")
       }
       else if (is.null(geneset) == FALSE){
         names(a)<-c("GE","Meth")
       }
       else if (is.null(snpset) == FALSE){
         names(a)<-c("SNP", "Meth")
       }
       else {
         names(a)<-c("Meth")
       }
     }
    
  invisible(capture.output(MOFAobject <- createMOFAobject(a)))
  ModelOptions <- getDefaultModelOptions(MOFAobject)
  DataOptions <- getDefaultDataOptions()
  # DirOptions <- list("dataDir" = tempdir(), "outFile" = tempfile()) ## old
  invisible(capture.output(TrainOptions <- getDefaultTrainOptions()))
  TrainOptions$DropFactorThreshold <- input$mofathresh
  TrainOptions$seed <- runif(1)
 # TrainOptions$maxiter <- 5
  invisible(capture.output(MOFAobject <- prepareMOFA(MOFAobject, DataOptions = DataOptions, ModelOptions = ModelOptions, TrainOptions = TrainOptions)))
  invisible(capture.output(MOFAobject <- try(runMOFA(MOFAobject))))
  #fil <- tempfile("/home/shiny/.virtualenvs/r-reticulate/mofaobject", fileext = ".rds")
  #saveRDS(MOFAobject, fil)
  #readRDS("/home/shiny/.virtualenvs/r-reticulate/mofaobject.rds")
  
  if ("try-error" %in% class(MOFAobject)) {     #(class(MOFAobject) == "try-error") {
    p <- 1
  }
  else{
   r2 <- calculateVarianceExplained(MOFAobject)
   # v factor reduction v
    # fact <- apply(r2$R2PerFactor[,c(names(a))], 1, sum)
    #for(m in 1:length(fact))
    #  if(fact[m]>=input$mofathresh){
    #    count <- m  
    #  }
    MOFAfactors <- getFactors(MOFAobject, factors="all", as.data.frame = F)
    dtcov <- data.frame(Group, MOFAfactors)
    if(input$select==2){
      model <- glm(Group~.,family=binomial(link='logit'),data=dtcov) #alt
      model2 <- glm(Group~1,family=binomial(link='logit'),data=dtcov) #null
    }
    else{
      model <- coxph(Surv(Observedtime, Group)~.,data=dtcov) 
      model2 <- coxph(Surv(Observedtime, Group)~1,data=dtcov)  
    }
    anova <- anova(model,model2,test="Chisq")  
    sumup <- anova
    p <- anova[2,5] 
  }

pdf[h,] <- p
 if(pdf[h,]<=input$error){
   value <- value + 1
 } 
   }
    #nmf + LR
     if(input$select==6){
       # SNF splitting as a list can also be used here
       GBM <- data.frame(matrix(ncol=0,nrow=n))
       if (is.null(geneset) == FALSE){
        GBM <- allGE
      }
      if (is.null(snpset)== FALSE){
        GBM <- cbind(GBM,allsnp)
      #  GBM <- append(GBM,list(allsnp))
      }
      if (is.null(methset)== FALSE){
        GBM <- cbind(GBM,allmeth)
      }
      if (input$tratio > 0){
        GBM <- cbind(GBM,Treatment)
      }
      GBM <- data.matrix(GBM)
     # GBM <- split(GBM, seq(ncol(GBM)))
     # GBM <- list(GBM)
      result=ExecuteCNMF(t(GBM),clusterNum=3,nrun=5)
      clust <- result$group
      model <- glm(Group~clust,family=binomial(link='logit'))
      sumup <- summary(model)
      p <- sumup
      pdf[h,] <- p$coefficients[2,4]
      if(pdf[h,]<=input$error){
        value <- value + 1
      } 
    }
    #snf + LR
    if(input$select==12){
      GBM <- c()
      if (is.null(geneset) == FALSE){
        GBM1 <- data.matrix(allGE)
        GBM <- list(t(GBM1))
      }
      if (is.null(snpset)== FALSE){
        GBM2 <- data.matrix(allsnp)
        GBM <- append(GBM,list(t(GBM2)))
      }
      if (is.null(methset)== FALSE){
        GBM3 <- data.matrix(allmeth)
        GBM <- append(GBM,list(t(GBM3)))
      }
     # if (input$tratio > 0){
     #  GBM <- cbind(GBM,Treatment)
     # }
      
      result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
      clust <- result$group
      model <- glm(Group~clust,family=binomial(link='logit'))
      sumup <- summary(model)
      p <- sumup
      pdf[h,] <- p$coefficients[2,4]
      if(pdf[h,]<=input$error){
        value <- value + 1
      } 
    }
    # survival snf
    if(input$select==13){
      GBM <- c()
      if (is.null(geneset) == FALSE){
        data1=FSbyCox(t(data.matrix(allGE)),Observedtime,Group,cutoff=0.05)
        GBM=list(data1)
      }
      if (is.null(snpset)== FALSE && is.null(geneset) == FALSE){
        data2=FSbyCox(t(data.matrix(allsnp)),Observedtime,Group,cutoff=0.05)
        GBM=append(GBM,list(data2))
      }
      else if (is.null(snpset)== FALSE && is.null(geneset) == TRUE) {
        data1=FSbyCox(t(data.matrix(allsnp)),Observedtime,Group,cutoff=0.05)
        GBM=list(data1)
      }
      if (is.null(methset)== FALSE && is.null(geneset) == FALSE && is.null(snpset) == FALSE){
        data3 <- FSbyCox(t(data.matrix(allmeth)),Observedtime,Group,cutoff=0.05)
        GBM=append(GBM,list(data3))
      }
      else if (is.null(methset)== FALSE && is.null(geneset) == TRUE && is.null(snpset) == TRUE) {
        data1=FSbyCox(t(data.matrix(allmeth)),Observedtime,Group,cutoff=0.05)
        GBM=list(data1)
      }
      else if ((is.null(methset)== FALSE && is.null(geneset) == TRUE && is.null(snpset) == FALSE) ||(is.null(methset)== FALSE && is.null(geneset) == FALSE && is.null(snpset) == TRUE) ){
        data2=FSbyCox(t(data.matrix(allsnp)),Observedtime,Group,cutoff=0.05)
        GBM=append(GBM,list(data2))
      }
     # if (input$tratio > 0){
    #  }
   
    result1<- try(ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20,plot = FALSE))
    if ("try-error" %in% class(result1)) {     #(class(MOFAobject) == "try-error") {
      pdf[h,] <- 1
    }
    else{
    group1=result1$group
    #LR <- survdiff(Surv(Observedtime, Group) ~ group1)
    coxmod <- coxph(Surv(Observedtime, Group) ~ group1)
    sumup <- summary(coxmod)
    p <- sumup
    pdf[h,] <- p$coefficients[1,5]
    }
    if(pdf[h,]<=input$error){
      value <- value + 1
    } 
    }
    
  # MCIA - correlation between omics p-value
    if(input$select==18){
      data <- c()
      if (is.null(geneset) == FALSE){
        dta <- list(t(allGE))
      }
      if (is.null(snpset)== FALSE){
        dta <- append(dta,list(t(allsnp)))
      }
      if (is.null(methset)== FALSE){
        dta <- append(dta,list(t(allmeth)))
      }
     
    mcoin <- mcia(dta)
    try(
    a <- cor.test(mcoin$mcoa$cov2[,1],mcoin$mcoa$cov2[,2])
    )
    sumup <- mcoin
    p <- a$p.value
    pdf[h,] <- p
    if(pdf[h,]<=input$error){
      value <- value + 1
    } 
    }
   pd[h,1]<- z*input$obs
    } #end of methods bracket - end of h
    # labels for plot
      colnames(pcaldf)[1] <- "SampleSize"
      colnames(pcaldf)[2] <- "Power"
      power <- value * 100 / (input$replicates)
      pcaldf[z,] <- list(n,power)#,row.names = "%")
      if(z == 1){
      pdf2 <- cbind(pd,pdf)
      }
      else{
      pdfadd <- cbind(pd,pdf)
      pdf2 <- rbind(pdf2,pdfadd)
      }
    }
    })
    melted = melt(pcaldf, id.vars="SampleSize")
    colnames(pdf2)[1] <- "SampleSize"
    colnames(pdf2)[2] <- "pvalue"
    
    #fdr page histogram
    output$distPlot <- renderPlotly({
      
     hist <-  ggplot(data = pdf2, aes(x=pvalue)) + geom_histogram() + xlab("p-value") + ylab("Frequency")
     ggplotly(hist) %>%
     layout(width = 800)#, legend = list(orientation = "h", x = 0.4, y = -0.2))
      # x    <- pdf2 #[,1]
      # 
      # hist(x, col = "#75AADB", border = "white",
      #      xlab = "p-values",
      #      main = "Histogram of p-values")
    })
    # FDR
    # Store p-values in a matrix
    exp.sig<- (input$replicates)*input$error
    obs.sig <- sum(pdf2[,2]<input$error)
  
    # Compare this with the Benjamini-Hochberg method:
    pvals.adj <- p.adjust(pdf2[,2], method="BH")
    dataset <- cbind(pdf2,pvals.adj)
    FDR <- exp.sig / obs.sig
    
    # Either 0.05 or bonferroni cutoff for # of genes
    output$false <- renderPrint({exp.sig / obs.sig})
    
    output$fdrplot<-renderPlotly({
    fdrgg <- ggplot(dataset,aes(x=pvalue,y=pvals.adj,color=as.factor(SampleSize)))+ geom_point(size=2)+labs(x = "p-values",y="Benjamini-Hochberg Adjusted p-values")+ geom_vline(xintercept = input$error,color="red") + labs(color = "")
    ggplotly(fdrgg) %>%
    layout(width = 800)#,legend = list(orientation = "h", x = 0.4, y = -0.2))
    })
    
    output$plot<-renderPlotly({
    gg <- ggplot(melted,aes(x=SampleSize,y=value, group=variable, colour = variable))+ geom_point(size=3)+geom_line(size=2)+labs(x = "Sample Size",y="Power %",colour="Omics Variable")
    gg2 <- gg + scale_colour_discrete(name="Omics Variable",labels=dput(names(pcaldf))) #+ scale_color_manual(values=c("forestgreen", "steelblue3"),labels=c("SNPs","Gene Expression"))
    ggplotly(gg2)
    
    })
    
    output$powersummary <- renderPrint({pcaldf})
    sumup
    })
   # variable selection check box group
  observe({
    x <- input$omics
   choices=c()
   if ('Genome' %in% x) 
    choices = c(choices,list("SNP-set" = 1))
    if("RNA-seq" %in% x) 
      choices =c(choices, list("GE" = 2))
    if("Epigenome" %in% x) 
      choices = c(choices,list("Meth" = 3))
    if(input$tratio > 0) 
      choices = c(choices,c("Treatment" = 4))
    
    updateCheckboxGroupInput(session, "checkGroup",
                             label = paste("Parameters"),
                             choices=choices,
                          #   selected = choices
                             )
    
    if(input$select == 5 || input$select == 17){
      if (is.null(choices))
        choices <- character(0)

      updateSelectInput(session, "mediator",
                        label = "Select mediator variable",
                        choices = choices,
                     #   selected = tail(choices, 1)
                     )
      updateSelectInput(session, "causalcov",
                        label = "Select dependent causal variable",
                        choices = choices,
                     #   selected = tail(choices, 2)
      )
    }
  })
  
  # Notifications
   observeEvent(6 %in% input$select, {
    if(6 %in% input$select){
     showNotification("NMF is very computationally demanding, therefore it may take longer to complete the calculation.",type = "warning")
     }
    })
   observeEvent(2 %in% input$select || 3 %in% input$select, {
     if(2 %in% input$select || 3 %in% input$select){
     showNotification("MOFA is very computationally demanding, therefore it may take longer to complete the calculation.",type = "warning")
     }
    })
   
  output$summary <- renderPrint({methods()})
  output$meanmeth <- renderText({paste0("Cases mean methylation ", round(mean(rbeta(input$obs,input$casemeth,input$casemeth2)),3),"\nControls mean methylation ", round(mean(rbeta(input$obs,input$contmeth,input$contmeth2)),3))
  })
 # output$powercalc <- renderText({ 
 #   paste0("In order to achieve ", input$target, "% power, you'll need to use a sample size of at least 50.")
 # })
  output$dloadData <- downloadHandler(
    filename = function() {
      paste("data-", ".txt", sep = "")
    },
    content = function(file) {
      write.table(datainputs(), file, row.names = FALSE)
    }
  )

})
