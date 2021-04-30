clr  = function()
{
  rm(list = ls(envir = .GlobalEnv),
     envir = .GlobalEnv)
}
clr()

library(shiny)
library(shinythemes)
library(readxl)
library(factoextra)
library(gridExtra)
library(GGally)
library(ggpubr)
library(cowplot)
library(SomaDataIO)
library(moments)
library(data.table)
library(DT)
source("shinyPlotSource.R")

### initialisation for all the user selected inputs, which are not reactive
inputfile1 <- "/Users/ydeng/Documents/QCstepOA/SS-200008.ADat"
inputfile2 <- "/Users/ydeng/Documents/QCstepOA/SS-200008.hybNorm.medNormInt.plateScale.medNormRefSMP.ADat"   
inputfile3 <- "/Users/ydeng/Documents/QCstepOA/STEpUP_QCData_Tranche1.xlsx"
inputfile4 <- "/Users/ydeng/Documents/QCstepOA/1_Master list_analyte concentrations.xlsx"

RawM <- initQCnorm(inputfile1,inputfile2) ###raw RFUs from Adat, ajust SampleType according to tranch Excel
RawMList = ExtractClinicG(RawM) ### adjust columns for our own case
RawM = RawMList[[1]]
MetaRawM = RawMList[[2]]

###for External validation
toTest1 <- c("mcp1bl","il6bl","il8bl","mmp3bl","activinabl","tsg6bl","timp1bl","tgfb1bl","fgf2bl")
toTest2 <- c("CCL2.2578.67","IL6.4673.13","CXCL8.3447.64","MMP3.2788.55","INHBA.13738.8","TNFAIP6.5036.50","TIMP1.2211.9","TGFB1.2333.72","FGF2.3025.50")
SandwichExtract <- ExtractSandwich(inputfile4)
sandwich_master_Ben = SandwichExtract[[1]] 
sandwich_master_Historic = SandwichExtract[[2]] 


###set up UI
shinyUI <- fluidPage(theme = shinytheme("cerulean"),
                     navbarPage("SomaLogic Proteomics Quality Control App",
                                tabPanel("Analysis Zone",
                                         titlePanel(h4("Proteomics Data Normalisation and Evaluation Methods")),
                                         sidebarLayout(position="left",
                                                       sidebarPanel(fileInput("uploadfile", "Upload SomaLogic Adat file"),
                                                                    textInput("norms","Select normalisation steps",""),
                                                                    radioButtons("Assessment","Select the assessment",
                                                                                 list("Total Protein","Limit of Detection","Calibrator Check","Coefficient Variance Breakdown","External Validation","KNN Batch Effect Test","Anova per Protein"),""),
                                                                    selectInput("PCAx", "Choose a parameter for PCA display:",
                                                                                list("","PlateID","Disease Group","Corhort","Blood Stain","Sample Age","Variation Explained by Principle Components"))),
                                                       
                                                       
                                                       mainPanel(h4("Assessment result display"),
                                                                 tabsetPanel(
                                                                   tabPanel("Plot for Radio Button", plotOutput("ptFinal1",width="100%",height="750px")), 
                                                                   tabPanel("Plot for Select Box", plotOutput("ptFinal2",width="100%",height="750px")), 
                                                                   tabPanel("Table", tableOutput("tblFinal"))),
                                                                 downloadButton("download1")
                                                                 )
                                                       )
                                         ),
                                tabPanel("User Guide")
                     )
)


###set up server
options(shiny.maxRequestSize=30*1024^2)
shinyServer <- function(input,output){
  
  ### define MySoma ect. as reactive variables, save computational cost
  
  ###output plots  
  output$ptFinal1 <- renderPlot({
    MySoma <- setReactiveMySoma(input$norms,RawM)
    AdatExtract <- ExtractAdat(MySoma,MetaRawM)
    exprDat_norm =  AdatExtract[[1]]
    calib_norm = AdatExtract[[2]]
    plasma_norm = AdatExtract[[3]]
    calibIDs = AdatExtract[[4]]
    calibPlates = AdatExtract[[5]]
    MetaRaw = AdatExtract[[6]]
    TranchExtract <- ExtractTranch(inputfile3,MySoma)
    metadata_reord = TranchExtract[[1]] 
    bloodStain = TranchExtract[[2]] 
    sampleAge = TranchExtract[[3]] 
    patientAge = TranchExtract[[4]] 
    plateID = TranchExtract[[5]]
    
    if(input$Assessment == "Total Protein"){TotalProCheck(exprDat_norm,bloodStain,metadata_reord)}
    
    if(input$Assessment == "Calibrator Check"){CalibratorCheckPlot(calib_norm,calibIDs)}
    
    if(input$Assessment =="Coefficient Variance Breakdown"){CVbreakPlot(calib_norm,calibIDs,calibPlates)}
    
    if(input$Assessment == "Limit of Detection"){LoDdetection(RawM)}
    
    if(input$Assessment=="Anova per Protein"){ConfounderTable <- ConfouderCheck(exprDat_norm,plateID,metadata_reord,bloodStain,sampleAge)
    confounderPlot2(ConfounderTable)}
  })
  
  output$ptFinal2 <- renderPlot({
    MySoma <- setReactiveMySoma(input$norms,RawM)
    AdatExtract <- ExtractAdat(MySoma,MetaRawM)
    exprDat_norm =  AdatExtract[[1]]
    calib_norm = AdatExtract[[2]]
    plasma_norm = AdatExtract[[3]]
    calibIDs = AdatExtract[[4]]
    calibPlates = AdatExtract[[5]]
    MetaRaw = AdatExtract[[6]]
    TranchExtract <- ExtractTranch(inputfile3,MySoma)
    metadata_reord = TranchExtract[[1]] 
    bloodStain = TranchExtract[[2]] 
    sampleAge = TranchExtract[[3]] 
    patientAge = TranchExtract[[4]] 
    plateID = TranchExtract[[5]]
    ptPCAfinal <- PCAglob(MetaRaw,exprDat_norm,input$PCAx)
    ptPCAfinal
  })
  
  ###output tables  
  output$tblFinal <- renderTable({
    MySoma <- setReactiveMySoma(input$norms,RawM)
    AdatExtract <- ExtractAdat(MySoma,MetaRawM)
    exprDat_norm =  AdatExtract[[1]]
    calib_norm = AdatExtract[[2]]
    plasma_norm = AdatExtract[[3]]
    calibIDs = AdatExtract[[4]]
    calibPlates = AdatExtract[[5]]
    MetaRaw = AdatExtract[[6]]
    TranchExtract <- ExtractTranch(inputfile3,MySoma)
    metadata_reord = TranchExtract[[1]] 
    bloodStain = TranchExtract[[2]] 
    sampleAge = TranchExtract[[3]] 
    patientAge = TranchExtract[[4]] 
    plateID = TranchExtract[[5]]
    
    if(input$Assessment == "External Validation"){
      CorData_norm1 <- ExtVal1(sandwich_master_Ben,exprDat_norm,toTest1,toTest2,metadata_reord,"OA")
      CorData_norm2 <- ExtVal1(sandwich_master_Historic,exprDat_norm,toTest1,toTest2,metadata_reord,"Knee Injury")
      CorData = as.table(as.matrix(rbind(CorData_norm1,CorData_norm2)))
      CorData}
    
    if(input$Assessment == "KNN Batch Effect Test"){BatchEffectM <- KNNtest(exprDat_norm,MetaRaw)
    BatchEffectM}
  },colnames=FALSE)
}

shinyApp(ui = shinyUI, server = shinyServer)
