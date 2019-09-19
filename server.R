# MOPower v1.0.2
# This is the server logic of the Shiny web application. You can run the 
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
library(MOFA)
library(edgeR)
library(foreach)
library(doParallel)

# v This section is for shinyapps.io v
#use_python("/usr/bin/python3.5")
#virtualenv_create("r-reticulate")
#use_virtualenv("r-reticulate")
#py_install("mofapy")

# Define server logic

#plan(multisession)

graph1 <- as.data.frame(matrix(NA, ncol = 0, nrow = 0))

shinyServer(function(input, output, session) {
  
  # Reactive datatable simulation
  # Simulation of download data and UI display data
  
 datainputs <- reactive({
  # Study design parameters
  # Reactive datatable
  n <- input$obs
  EndofStudy <- input$eos
  Treatment <- rbinom(n,1,input$tratio)
  # Defining gene-level parameters
  numsnp <- input$snps
  numgene <- input$genes
  Patients <- 1:n
  sum1 <- -0
  dt <- data.frame(Patients,Treatment)
  showdt <- data.frame(Patients,Treatment)
  
  # SNPs and rare variants
  if("Genome" %in% input$omics){
    snpdt <- as.data.frame(matrix(NA, ncol = numsnp, nrow = n))
    beta <- c(1:numsnp)
    agg<- as.data.frame(matrix(NA, ncol = 4, nrow = n))
    linpred <- matrix(NA, nrow = n,ncol=numsnp)
    # SNP data simulation
    for(i in 1:numsnp){
      maf <- runif(1,input$mafint[1], input$mafint[2])
      snpdt[i] <- cbind(rbinom(n,2,maf))
      colnames(snpdt)[i] <- paste("SNP", i, sep = "")
      beta[i] <- sqrt((runif(1,input$snpeffect[1], input$snpeffect[2]))/2/maf/(1-maf))
    }
    agg <- t(t(snpdt)*beta)
    sum1 <- apply(agg, 1, sum)
    dt <- cbind(dt,snpdt)
    showdt$SNP <- snpdt[,1]
  }
  
  # Outcome (case-control)
  if(1 %in% input$outselect){
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
    showdt$Group <- Group
  }
  
  # Outcome (Time-to-event)
  else if(3 %in% input$outselect){
    # baseline hazard
    lambda <- input$Stime*exp(-sum1 - (Treatment*input$treat))
    survtime <- rweibull(input$obs,1, scale=lambda)
    censtime <- rbinom(input$obs,1,1-(input$Ctime*0.01))
    # Right censoring as percentage
    Group <- ifelse(survtime >= input$eos,0,censtime)
    Observedtime <- ifelse(survtime>=input$eos,input$eos,survtime)
    # v Censoring is random within sample v
    #Observedtime <- ifelse(survtime>input$eos & censtime>input$eos,input$eos,pmin(survtime,censtime))
    #Group <- ifelse(survtime < censtime & survtime < input$eos,1,0)
    dt$Censoring <- Group
    dt$Observedtime <- Observedtime
    dt$EndofStudy <- EndofStudy
    showdt$Censoring <- Group
    showdt$Observedtime <- Observedtime
    showdt$EndofStudy <- EndofStudy
  }
  
  # Outcome (longitudinal) # needs to be a new gene expression value for each timepoint
  #if(2 %in% input$outselect){
    #timevec <- matrix(rep(seq(0, input$eos, input$fusamples), input$obs), ncol=(input$eos/input$fusamples)+1, nrow=input$obs, byrow=T)
    #colnames(timevec) <- colnames(timevec, do.NULL = FALSE, prefix = "Time")
    #dt <- cbind(dt,timevec)
    #showdt$TimePoint <- timevec[,1]
  #}
  
  # Gene expression data simulation
  if("RNA-seq" %in% input$omics){
    mexp <-  exp(Group*(runif(1,input$fc[1],input$fc[2]))+log(input$reads)) # N can be input or estimated from uniform distribution. N = library size
    disp <- 1/(runif(1,input$disp[1],input$disp[2]))#^2 #  or variance = mexp*(1+mexp*input$disp)or variance=mexp+(mexp^2/input$disp)
    GE <- as.data.frame(matrix(NA, ncol = numgene, nrow = n))
    for(j in 1:numgene){
    GE[j] <- cbind(rnbinom(n, size=disp, mu=mexp))
    colnames(GE)[j] <- paste("GE", j, sep = "")
    }
  # DT without Methylation
  dt <- cbind(dt,GE)
  showdt$GeneExp <- GE[,1]
  }
  
  # DNA methylation - array CpG sites
  if("Epigenome" %in% input$omics){
    sites <- input$cpg
    # r and v are unmeth and hemimeth column numbers
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
        else {
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
    showdt$CpGMeth <- methsites[,1]
  }
  
  # Sample data and download data need to be seperated in a list
  list(a=dt,b=showdt)
 })
  
  # Render data in UI
 output$table <- renderTable({
  head(datainputs()$b, n = 20)
 })
 
##################################### End of sample data and download data simulation ############################################# 
 
 # Power calculation section
 # Simulate data, integration analysis models and power
 #methods <- eventReactive(input$submit, {
 methods <- observeEvent(input$submit, {
  
  # Dataframes for storing p-values and sample size number
  pcaldf <- as.data.frame(matrix(NA, ncol = 2, nrow = 10))
  pdf2 <- as.data.frame(matrix(NA, ncol = 2, nrow = 1))
  pdfadd <- as.data.frame(matrix(NA, ncol = 2, nrow = 1))
  estdisp <- as.data.frame(matrix(NA, ncol = 1, nrow = 0))
  
  # Number of simulations
  simval <- ifelse("RNA-seq" %in% input$omics && ((19 %in% input$select) || (20 %in% input$select)),1,input$replicates)

  # All conditions outside loop
  cond1 <- "Genome" %in% input$omics
  cond2 <- 1 %in% input$outselect
  cond3 <- 3 %in% input$outselect
  cond4 <- "RNA-seq" %in% input$omics && ((19 %in% input$select) || (20 %in% input$select))
  cond5 <- "Epigenome" %in% input$omics
  cond6 <- 1 %in% input$checkmeth
  cond7 <- 2 %in% input$checkmeth
  cond8 <- input$select==1
  cond12 <- input$tratio > 0 && 4 %in% input$checkGroup
  cond13 <- input$select==16
  cond14 <- input$select==19
  cond15 <- input$select==20
  cond16 <- input$select==10
  cond17 <- input$select==5 || input$select == 17
  cond18 <- input$select==2 || input$select==3
  cond19 <- input$select==6
  cond20 <- input$select==12
  cond21 <- input$select==13
  cond22 <- input$select==18
  
  # shiny variables
  nn <- input$obs
  ee <- input$eos
  tt <- input$tratio
  ns <- input$snps
  ng <- input$genes
  mf1 <- input$mafint[1]
  mf2 <- input$mafint[2]
  snpe1 <- input$snpeffect[1]
  snpe2 <- input$snpeffect[2]
  cc <- input$ccratio
  trt <- input$treat
  stt <- input$Stime
  ctt <- input$Ctime
  fcc <- input$fc
  rds <- input$reads
  dss <- input$disp
  unm <- input$unmeth
  hmm <- input$hemimeth
  cpg1 <- input$cpg
  cmm <- input$casemeth
  cmm2 <-input$casemeth2
  ctm <- input$contmeth
  ctm2 <- input$contmeth2
  err <- input$error
  sel <- input$select
  mthresh <- input$mofathresh
  med <- input$mediator
  ccov <- input$causalcov
  check <- input$checkGroup
  # Progress bar
   withProgress(message = 'Calculating Power', value = 0, {
  
  cl <- makeCluster(5)
  registerDoParallel(cl) 
     
  # Power for 10 different sample sizes (User specified sample size --> increasing)
    for (z in (1:10)){ 
    #  future({
      n <- nn*z
      value <- 0
      value2 <- 0
      # Matrix to store p-values
      pdf <- as.data.frame(matrix(NA, ncol = 1, nrow = 1))
      # Matrix to store sample size number
      pd <- as.data.frame(matrix(NA, ncol = 1, nrow = 0))
      new <- as.data.frame(matrix(NA, ncol = 3, nrow = simval))
      
      # Progress bar
      incProgress(1/z, detail = paste("Please wait...."))

      # Simulate all data replicates  
      paral <- foreach(h=1:simval, .combine=rbind, .packages = c("shiny","shinyBS","edgeR","MOFA","survival","SNFtool", "iCluster", "iClusterPlus","omicade4","NMF","mediation","ConsensusClusterPlus","glmpath","lme4","extRemes","CancerSubtypes")) %dopar% {
     #   for (h in 1:simval){
       n <- nn*z
        sites <- 0
        unsites <- 0
        hemisites <- 0
        EndofStudy <- ee
        methset <- NULL
        geneset <- NULL
        snpset <- NULL
        sum1 <- 0
        Treatment <- rbinom(n,1,tt)
        numsnp <- ns
        numgene <- ng
        Patients <- 1:n
        clnum <- 0
        
        # SNPs and rare varaints
        if(cond1){
          clnum <- clnum + 1
          snpdt <- as.data.frame(matrix(NA, ncol = numsnp, nrow = n))
          beta <- c(1:numsnp)
          agg<- as.data.frame(matrix(NA, ncol = 4, nrow = n))
          linpred <- matrix(NA, nrow = n,ncol=numsnp)
          for(i in 1:numsnp){
            maf <- runif(1,mf1, mf2)
            snpdt[i] <- cbind(rbinom(n,2,maf))
            colnames(snpdt)[i] <- paste("SNP", i, sep = "")
            beta[i] <- sqrt((runif(1,snpe1, snpe2)/2/maf/(1-maf)))
          }
          agg <- t(t(snpdt)*beta)
          sum1 <- apply(agg, 1, sum)
          dt <- data.frame(Patients,Treatment,EndofStudy,snpdt)
          allsnp <- dt[ , grepl( "SNP" , names( dt ) ) ]
          snpset <- apply(allsnp, 1, sum)
        }
        else{
          dt <- data.frame(Patients,Treatment,EndofStudy)
        }
      
        # Outcome (case-control)
        if(cond2){
          if(cond1){
            b0 <- log(cc/(1-cc)) -  sum1 - (Treatment*trt)
            linpred <- b0 + sum1 + (Treatment*trt)
            pi <- exp(linpred) / (1 + exp(linpred))
            Group <- rbinom(n, size=1, prob=pi)
          }
          else{
            Group <- rbinom(n, size=1, prob=cc)
          }
          dt$Group <- Group
        }
      
        # Outcome (survival)
        else if(cond3){
          # baseline hazard
          lambda <- stt*exp(-sum1 - (Treatment*trt))
          survtime <- rweibull(nn*z,1, scale=lambda)
          censtime <- rbinom(n,1,1-(ctt*0.01))
          # right censoring
          Group <- ifelse(survtime >= ee,0,censtime)
          Observedtime <- ifelse(survtime>=ee,ee,survtime)
          dt$Group <- Group
          dt$Observedtime <- Observedtime
        }
      
        # Gene expression
        # For exact test and NB model, data is simulated for multiple genes for each individual.
        if(cond4){
          value<-0
          clnum <- clnum + 1
          GE <- as.data.frame(matrix(NA, ncol = n, nrow = numgene))
          for (i in 1:n){
            if (Group[i] == 1){
              mexp <-  exp(1*fcc+log(rds))
            }
            else{
              mexp <-  exp(0*fcc+log(rds))
            }
            GE[i]<-cbind(rnbinom(numgene, size=1/(dss), mu=mexp))
          }
          # edgeR format conversion
          y <- DGEList(counts=GE, group=Group)
          y <- calcNormFactors(y)
          design <- model.matrix(~Group)
          y <- estimateDisp(y,design)
          di <- estimateDisp(y, design)
          d1 <- as.data.frame(di$tagwise.dispersion) 
          estdisp <- rbind(estdisp,d1)
        }
        # Single or multiple gene simulation. Transpose of above simulation.
        else{
          clnum <- clnum + 1
          mexp <-  exp(Group*fcc+log(rds)) 
          disp <- 1/(dss) 
          GE <- as.data.frame(matrix(NA, ncol = numgene, nrow = n))
          for(j in 1:numgene){
            GE[j] <- cbind(rnbinom(n, size=disp, mu=mexp))
            colnames(GE)[j] <- paste("Genecount", j, sep = "")
          }
          dt <- cbind(dt,GE)
          allGE <- dt[ , grepl( "Genecount" , names( dt ) ) ]
          geneset <- apply(allGE, 1, sum)
        }
     
        # DNA methylation - array CpG sites
        if(cond5){
          sites <- cpg1
          clnum <- clnum + 1
          v <- 0
          r <- 0
          unsites <- as.data.frame(matrix(NA, ncol = unm, nrow = n))
          hemisites <-as.data.frame(matrix(NA, ncol = hmm, nrow = n))
          methsites <- as.data.frame(matrix(NA, ncol = sites, nrow = n))
          for(k in 1:sites){
            meth = ifelse(Group==1,rbeta(1,cmm,cmm2),rbeta(1,ctm,ctm2))
            methsites[,k] = cbind(meth)
            colnames(methsites)[k] <- paste("CpGmeth", k, sep = "")
          }
          # unmeth and hemi meth sites
          if(cond6){
            for(v in 1:unm){
            unsites[,v] <- rbeta(n,1,6)
            colnames(unsites)[v] <- paste("CpGmeth", v+k+r, sep = "")
            }
            methsites <- cbind(methsites,unsites)
          }
          if(cond7){
            for(r in 1:hmm){
              hemisites[,r] <- rbeta(n,6,6)
              colnames(hemisites)[r] <- paste("CpGmeth", r+v+k, sep = "")
            }
            methsites <- cbind(methsites,hemisites)
          }
          dt <- cbind(dt,methsites)
          allmeth <- dt[ , grepl( "CpGmeth" , names( dt ) ) ]
          methset <- apply(allmeth, 1, sum)
        }
      
        
        # Analysis methods section
        
        # Logistic regression model
        # Data formatting
        if(cond8){
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
          if (cond12){
            dta <- cbind(dta,Treatment)
          }
        
          # GLM model
          model <- glm(Group~.,data=dta,family=binomial(link='logit'))
          # sumup is printed in UI as an analysis replicate
          #sumup <- summary(model)
          model2 <- glm(Group~1,data=dta,family=binomial(link='logit'))
          anova <- anova(model,model2,test="Chisq")  
          pval1 <- anova[2,5]
          pdf[h,] <- pval1

        }
  
        # Cox peoportional hazards model  
        if(cond13){
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
          if (cond12){
            dta <- cbind(dta,Treatment)
          }
      
          model <- coxph(Surv(Observedtime, Group) ~.,data=dta)
          #sumup <- summary(model)
          model2 <- coxph(Surv(Observedtime, Group) ~1,data=dta)
          anova <- anova(model,model2,test="Chisq")  
          pval1 <- anova[2,5]
          pdf[h,] <- pval1

        }
    
        # Exact test (RNA-seq data only)
        if(cond14){
          et <- exactTest(y)
          #sumup <- topTags(et)
          for (k in 1:numgene){
            pdf[k,] <- et$table[k,3]
            if(pdf[k,]<=err){
              value <- value + 1
            }
          pd[nrow(pd)+1,1]<- z*nn
          }
          }
        
    
        # NB GLM
        if(cond15){
      fit <- glmQLFit(y, design)
      nbqlf <- glmQLFTest(fit,coef=2)
      #sumup <- topTags(nbqlf)
      for (k in 1:numgene){
        pdf[k,] <- nbqlf$table[k,4]
        if(pdf[k,]<=err){
          value <- value + 1
        }
       pd[nrow(pd)+1,1]<- z*nn
      }
    }
    
        # Cox Path L1 Penalty
        if(cond16){
      
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
      if (cond12){
        features <- cbind(features,Treatment)
      }
      dat <- list(x=features, time=Observedtime, status=Group)
      invisible(capture.output(fit.a <- coxpath(dat)))
      len <- length(fit.a$loglik)
      last <- tail(fit.a$df, n=1)
      b <- lr.test(fit.a$loglik[len], fit.a$loglik[1], alpha = err, df = last)
      pdf[h,] <- b$p.value
      #if(z==10 && h==simval){
      # sumup <- summary(fit.a) 
      #}
    }
    
      # Mediation analysis
      if(cond17){
      
      if(2 %in% med){
        mediator <- geneset
      }
      else if(1 %in% med){
        mediator <- snpset
      }
      else {
        mediator <- methset
      }
      
      if(2 %in% ccov){
        abc <- geneset
      }
      else if(1 %in% ccov){
        abc <- snpset
      }
      else {
        abc <- methset
      }
      
      dta <- data.frame(Group,mediator,abc)
      
      if(1 %in% check && ((2 %in% med && 3 %in% ccov) || (3 %in% med && 2 %in% ccov))){
        b <- snpset
        dta <- cbind(dta,b)
        }
      else if(2 %in% check && ((1 %in% med && 3 %in% ccov)|| (3 %in% med && 1 %in% ccov))){
        b <- geneset
        dta <- cbind(dta,b)
      }
      else if(3 %in% check && ((2 %in% med && 1 %in% ccov) || (1 %in% med && 2 %in% ccov))){
        b <- methset
        dta <- cbind(dta,b)
      }

      if (4 %in% check){
      c <- Treatment
      dta <- cbind(dta,c)
      }
      
      med.fit <- lm(mediator ~.,data=dta[,!(names(dta) %in% c('Group'))])
      if(sel==5){
      out.fit <- glm(Group~.*., data=dta,family = binomial("probit"))
      }
      else{
      out.fit <- survreg(Surv(Observedtime,Group)~.*., data = dta, dist='weibull')
      }
      med.out <- try(mediate(med.fit, out.fit, treat = "abc" , mediator = "mediator",robustSE = TRUE, sims = 100))
      if ("try-error" %in% class(med.out)) {
        pdf[h,] <- 1
      }
      else{
      #sumup <- summary(med.out)
      a <- summary(med.out)
      # total effect p-value
      pdf[h,] <- a$tau.p
      }
      
    }

    # MOFA
    if(cond18){
  
     a <- list()
     
     if (is.null(geneset) == FALSE){
       if(h==1 && z==1){
         rownames(allGE) <- paste("Sample", 1:nn*z, sep = "")
       }
       a <- list(t(allGE))
       names(a)<-c("GE")
     }
     if (is.null(snpset)== FALSE){
       if(h==1 && z==1){
         rownames(allsnp) <- paste("Sample", 1:nn*z, sep = "")
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
         rownames(allmeth) <- paste("Sample", 1:nn*z, sep = "")
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
  invisible(capture.output(TrainOptions <- getDefaultTrainOptions()))
  TrainOptions$DropFactorThreshold <- mthresh
  TrainOptions$seed <- runif(1)
  invisible(capture.output(MOFAobject <- prepareMOFA(MOFAobject, DataOptions = DataOptions, ModelOptions = ModelOptions, TrainOptions = TrainOptions)))
  invisible(capture.output(MOFAobject <- try(runMOFA(MOFAobject))))

  if ("try-error" %in% class(MOFAobject)) {     #(class(MOFAobject) == "try-error") {
    p <- 1
  }
  else{
   r2 <- calculateVarianceExplained(MOFAobject)

    MOFAfactors <- getFactors(MOFAobject, factors="all", as.data.frame = F)
    dtcov <- data.frame(Group, MOFAfactors)
    if(sel==2){
      model <- glm(Group~.,family=binomial(link='logit'),data=dtcov) #alt
      model2 <- glm(Group~1,family=binomial(link='logit'),data=dtcov) #null
    }
    else{
      model <- coxph(Surv(Observedtime, Group)~.,data=dtcov) 
      model2 <- coxph(Surv(Observedtime, Group)~1,data=dtcov)  
    }
    anova <- anova(model,model2,test="Chisq")  
    #sumup <- anova
    p <- anova[2,5] 
  }

pdf[h,] <- p

rm(MOFAobject)

   }
  
       # NMF + LR
        if(cond19){
       # SNF splitting as a list can also be used here
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
          
        #  GBM <- data.frame(matrix(ncol=0,nrow=n))
       #if (is.null(geneset) == FALSE){
      #  GBM <- allGE
     # }
     # if (is.null(snpset)== FALSE){
       # GBM <- cbind(GBM,allsnp)
      #  GBM <- append(GBM,list(allsnp))
    #  }
     # if (is.null(methset)== FALSE){
     #   GBM <- cbind(GBM,allmeth)
     # }
     # if (tt > 0){
       # GBM <- cbind(GBM,Treatment)
     # }
     # GBM <- data.matrix(GBM)
     # GBM <- split(GBM, seq(ncol(GBM)))
     # GBM <- list(GBM)
      result <- try(ExecuteCNMF(GBM,clusterNum=clnum,nrun=5))
      if ("try-error" %in% class(result)) {
        pdf[h,] <- 1
      }
      else{
      clust <- result$group
      model <- glm(Group~clust,family=binomial(link='logit'))
      p <- summary(model)
      pdf[h,] <- p$coefficients[2,4]
}
    }
  
        # SNF + LR
        if(cond20){
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
     # if (tt > 0){
     #  GBM <- cbind(GBM,Treatment)
     # }
     # par(mar = rep(2, 4))
      result <- try(ExecuteSNF(GBM, clusterNum=clnum, K=20, alpha=0.5, t=20))
      if ("try-error" %in% class(result)) {
        pdf[h,] <- 1
      }
      else{
      clust <- result$group
      model <- glm(Group~clust,family=binomial(link='logit'))
      p <- summary(model)
      pdf[h,] <- p$coefficients[2,4]
}
    }
    
        # SNF + Cox
        if(cond21){
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
   
    result1<- try(ExecuteSNF(GBM, clusterNum=clnum, K=20, alpha=0.5, t=20,plot = FALSE))
    if ("try-error" %in% class(result1)) {     #(class(MOFAobject) == "try-error") {
      pdf[h,] <- 1
    }
    else{
    group1=result1$group
    #LR <- survdiff(Surv(Observedtime, Group) ~ group1)
    coxmod <- coxph(Surv(Observedtime, Group) ~ group1)
    p <- summary(coxmod)
    pdf[h,] <- p$coefficients[1,5]
    }
    }
    
        # MCIA - correlation between omics p-values
      #  if(cond22){
     # dta <- c()
      #if (is.null(geneset) == FALSE){
      #  dta <- list(t(allGE))
      #}
      #if (is.null(snpset)== FALSE){
      #  dta <- append(dta,list(t(allsnp)))
      #}
      #if (is.null(methset)== FALSE){
      #  dta <- append(dta,list(t(allmeth)))
      #}
     #mcoin <- try(mcia(dta))
     #if ("try-error" %in% class(mcoin)) {
     #   pdf[h,] <- 1
      #}
      #else{
   # a <- try(cor.test(mcoin$mcoa$cov2[,1],mcoin$mcoa$cov2[,2]))
    #if ("try-error" %in% class(a)) {
    #  pdf[h,] <- 1
   # }
    #else{
         # sumup <- mcoin
   # p <- a$p.value
    #pdf[h,] <- p
   # }
    

      #}
  #  }
  
        # Populating sample size dataframe
        if(sel != 19 && sel != 20){
          pd[h,1]<- z*nn 
        }
        
        new <- data.frame(a=pdf[h,],b=pd[h,1])
      
        }
      
####################################### End of number of simulations ######################      
      
     
       value <- ifelse(paral[,1]<=err,value + 1,value)
      value <- sum(value, na.rm=T)
      
      # for (h in 1:simval) {
      #  if(paral[h,1]<=err){
      #  value <- value + 1
      # }
      #}

      # Labels for power text output
      colnames(pcaldf)[1] <- "SampleSize"
      colnames(pcaldf)[2] <- "Power"
      
      # number of sig genes or simulations
      power <- ifelse("RNA-seq" %in% input$omics && ((19 %in% input$select) || (20 %in% input$select)),value*100/(input$genes), value * 100 / (simval))
      pcaldf[z,] <- list(n,power)
      if(z == 1){
        pdf2 <- cbind(paral[,2],paral[,1])
      }
      else{
        pdfadd <- cbind(paral[,2],paral[,1])
        pdf2 <- rbind(pdf2,pdfadd)
      }
   #   })
      }
  
################### End of different sample size power calculation ###########################
  
  
  
    # Labels for p-value distribution plot
    melted = melt(pcaldf, id.vars="SampleSize")
    pdf2 <- data.frame(pdf2)
    colnames(pdf2)[1] <- "SampleSize"
    colnames(pdf2)[2] <- "pvalue"
    
    # p-value histogram
    output$distPlot <- renderPlotly({
      hist <-  ggplot(data = pdf2, aes(x=pvalue)) + geom_histogram(color="darkblue", fill="steelblue") + xlab("p-value") + ylab("Frequency")
      ggplotly(hist)
    })
    
    # FDR plot
    exp.sig<- (simval)*input$error
    obs.sig <- sum(pdf2[,2]<input$error)
  
    # Compare this with the Benjamini-Hochberg method:
    pvals.adj <- p.adjust(pdf2[,2], method="BH")
    dataset <- cbind(pdf2,pvals.adj)
    FDR <- exp.sig / obs.sig
    output$false <- renderPrint({exp.sig / obs.sig})
    output$fdrplot<-renderPlotly({
    fdrgg <- ggplot(dataset,aes(x=pvalue,y=pvals.adj,color=as.factor(SampleSize)))+ geom_point(size=2)+labs(x = "p-values",y="Benjamini-Hochberg Adjusted p-values")+ geom_vline(xintercept = input$error,color="red") + labs(color = "")# + theme(plot.margin = margin(2,.8,2,.8, "cm"))
    ggplotly(fdrgg)
    })
    
    # Dispersion plot
    if(cond4){
    colnames(estdisp)[1] <- "EstDispersion"
    output$disphist <- renderPlotly({
      hist1 <-  ggplot(data = estdisp, aes(x=EstDispersion)) + geom_histogram(color="forestgreen", fill="chartreuse3") + xlab("Estimated-Dispersion") + ylab("Frequency")
      ggplotly(hist1)
    })
    }
    
    output$powersummary <- renderPrint({pcaldf})
   # output$summary <- renderPrint({sumup})
    
    # Power plot - refresh and add line to plot
      if(input$addline==FALSE){
        graph1 <<- as.data.frame(matrix(NA, ncol = 0, nrow = 0)) 
      }
      var1 <- rep(input$label,10)
      graph <- melted
      graph <- cbind(graph,var1)
      graph1 <<- rbind(graph1,graph)
      
      if(input$addline==TRUE){
        output$plot<-renderPlotly({
          
          gg <- ggplot(graph1,aes(x=SampleSize,y=value, group=var1, colour = var1))+ geom_point(size=3)+geom_line(size=2)+labs(x = "Sample Size",y="Power %",colour=input$label)
          gg2 <- gg + scale_colour_discrete(name="",labels=input$label)
          ggplotly(gg2)
        })
      }
      else{
        output$plot<-renderPlotly({
          gg <- ggplot(graph,aes(x=SampleSize,y=value, group=var1, colour = var1))+ geom_point(size=3)+geom_line(size=2)+labs(x = "Sample Size",y="Power %",colour=input$label)
          gg2 <- gg + scale_colour_discrete(name="",labels=input$label)
          ggplotly(gg2)
        })
      }
      
  #    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"
    #    filename = "report.docx",
   #     content = function(file) {
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          #tempReport <- file.path(tempdir(), "report.rmd")
          #file.copy("report.rmd", tempReport, overwrite = TRUE)
          
          # Set up parameters to pass to Rmd document
          #params <- list(bb=pcaldf,a=input$snps)
            #a=input$obs, b=input$eos, c=input$tratio, d=input$snps,e=input$genes,f=input$mafint[1],
             #            g=input$mafint[2],h=input$snpeffect[1],i=input$snpeffect[2],j=input$ccratio,k=input$treat,l=input$Stime,
              #           m=input$Ctime,n=input$fc,o=input$reads,p=input$disp,q=input$unmeth,r=input$hemimeth,s=input$cpg,
               #          t=input$casemeth,u=input$casemeth2,v=input$contmeth,w=input$contmeth2,x=input$error,
                #         y=input$select,z=melted,aa=pdf2,bb=pcaldf)
          
          # Knit the document, passing in the `params` list, and eval it in a
          # child of the global environment (this isolates the code in the document
          # from the code in this app).
          #rmarkdown::render(tempReport, output_file = file,
         #                   params = params,
        #                    envir = new.env(parent = globalenv())
       #   )
    #    }
    #  )
    })
  })

############################## End of progress bar and methods section ###############################
  
 
  # Methods in implementation development 
  #"Generalised Linear Mixed Effects Model" = 4, 
  #"Clustering + LR" = 14, "Clustering + Cox" = 15, "Sparse partial least square regression model" = 7,
  #"Joint hierarchical Bayesian modelling (lme)" = 8, "GLM Path" = 9,
  #"Partial least squares (Mixomics)" = 11, "Multi Co-inertia analysis" = 18
  
  # Variable selection check box group & method update
  observe({
   x <- input$omics
   choices=c()
   choices2=c()
   
   if ('Genome' %in% x) 
    choices = c(choices,list("SNP-set" = 1))
    if("RNA-seq" %in% x)
      choices =c(choices, list("Gene Expression" = 2))
    if("Epigenome" %in% x) 
      choices = c(choices,list("Methylation" = 3))
    if(input$tratio > 0) 
      choices = c(choices,c("Treatment" = 4))
    
   updateCheckboxGroupInput(session, "checkGroup",
                             label = paste("Parameters"),
                             choices=choices
                             )
   
   # Mediation analysis selection choices
    if(input$select == 5 || input$select == 17){
      if (is.null(choices))
        choices <- character(0)

      updateSelectInput(session, "mediator",
                        label = "Select mediator variable",
                        choices = choices
                     )
      updateSelectInput(session, "causalcov",
                        label = "Select dependent causal variable",
                        choices = choices
      )
    }
  })
  
  # List of methods shown in UI, dependent on outcome selected
  observe({
    x <- input$omics
    choices2=c()
    if("RNA-seq" %in% x) {
      choices2 = c(choices2,list("Exact test" = 19, "Negative-Binomial GLM" = 20))
    }
    if(1 %in% input$outselect)
    choices2 = c(choices2,list("Joint regression model (LR)" = 1,"MOFA + Logistic regression" = 2,"Mediation analysis + LR" = 5,"Similarity Network Fusion + LR" = 12,"Joint Non-negative Matrix Factorization (iCluster+)" = 6))
    
    if(3 %in% input$outselect)
    choices2 = c(choices2,list("Joint regression model (Cox)"=16,"MOFA + Cox proportional hazards model" = 3,"Mediation analysis + Cox" = 17,"Cox Path L1 Penalty" = 10,"Similarity Network Fusion + Cox" = 13,"Joint Non-negative Matrix Factorization (iCluster+)" = 6))   
    
    if (is.null(choices2)){
      choices2 <- character(0)
    }
    else{
      updateSelectInput(session, "select", label = paste("Statistical integration/analysis models"),
                        choices=choices2)
    }
    })
  
  # Bonferroni correction value for the number of features
  output$BC <- renderPrint(
    {
      input$error2/input$genes2
    }
  )
  
  # Notifications
   observeEvent(6 %in% input$select, {
    if(6 %in% input$select){
     showNotification("NMF is very computationally demanding, therefore it may take longer to complete the calculation.",type = "warning",closeButton = TRUE)
     }
    })
   observeEvent(2 %in% input$select || 3 %in% input$select, {
     if(2 %in% input$select || 3 %in% input$select){
     showNotification("MOFA is very computationally demanding, therefore it may take longer to complete the calculation.",type = "warning",closeButton = TRUE)
     }
    })
   observeEvent("RNA-seq" %in% input$omics, {
     if("RNA-seq" %in% input$omics){
       showNotification("If performing single omics power calculation, simulate a single SNP, Gene or 
                        CpG site to analyse across multiple replicates. If performing a single omics 
                        power caculation of RNA-seq data, simulate multiple genes for Exact test and 
                        Negative Binomial GLM to achieve power based on the number of significantly 
                        associated differentially expressed genes. For other analyses input number 
                        of simulations and a single gene.",type = "warning",duration = 40, closeButton = TRUE)
     }
   })
   
   
   output$meanmeth <- renderText({paste0("Cases mean methylation ", round(mean(rbeta(input$obs,input$casemeth,input$casemeth2)),3),"\nControls mean methylation ", round(mean(rbeta(input$obs,input$contmeth,input$contmeth2)),3))
  })
  
  # Download data to file
   output$dloadData <- downloadHandler(
    filename = function() {
      paste("data", ".txt", sep = "")
    },
    content = function(file) {
      write.table(datainputs()$a, file, row.names = FALSE)
    }
  )
})
