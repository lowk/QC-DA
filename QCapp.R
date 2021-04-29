
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
source("shinyPlotSource.R")

### initialisation for all the user selected inputs, which are not reactive
inputfile1 <- "/Users/ydeng/Documents/QCstepOA/SS-200008.ADat"
inputfile2 <- "/Users/ydeng/Documents/QCstepOA/SS-200008.hybNorm.medNormInt.plateScale.medNormRefSMP.ADat"   
inputfile3 <- "/Users/ydeng/Documents/QCstepOA/STEpUP_QCData_Tranche1.xlsx"
inputfile4 <- "/Users/ydeng/Documents/QCstepOA/1_Master list_analyte concentrations.xlsx"

RawM <- initQCnorm(inputfile1,inputfile2) ###raw RFUs from Adat, adjust SampleType according to tranch Excel
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
                                                                   tabPanel("Plot for Radio Button", plotOutput("ptFinal1",width="100%",height="750px"),downloadButton("downloadPlot1","Download Plot")), 
                                                                   tabPanel("Plot for Select Box", plotOutput("ptFinal2",width="100%",height="750px"),downloadButton("downloadPlot2","Download Plot")), 
                                                                   tabPanel("Table", tableOutput("tblFinal"),downloadButton("downloadTable","Download Table")))
                                                       )
                                         )
                                ),
                                tabPanel("User Guide", h4("Data Group Lead by Luke Jostins"),
                                                          h5("Normalisation Steps Function Name: HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp,CALIBRATION"))
                     )
)

###set up server
options(shiny.maxRequestSize=30*1024^2)
shinyServer <- function(input,output){
  
  ### define MySoma ect. as reactive variables, save computational cost
  MySoma <-reactive(setReactiveMySoma(input$norms,RawM))
  AdatExtract <- reactive(ExtractAdat(MySoma(),MetaRawM))
  
  exprDat_norm =  reactive(AdatExtract()[[1]])
  calib_norm = reactive(AdatExtract()[[2]])
  plasma_norm = reactive(AdatExtract()[[3]])
  calibIDs = reactive(AdatExtract()[[4]])
  calibPlates = reactive(AdatExtract()[[5]])
  MetaRaw = reactive(AdatExtract()[[6]])
  
  TranchExtract <- reactive(ExtractTranch(inputfile3,MySoma()))
  metadata_reord = reactive(TranchExtract()[[1]])
  bloodStain = reactive(TranchExtract()[[2]])
  sampleAge = reactive(TranchExtract()[[3]])
  patientAge = reactive(TranchExtract()[[4]])
  plateID = reactive(TranchExtract()[[5]])
  
  
  ###output plots  
  output$ptFinal1 <- renderPlot({
    
    if(input$Assessment == "Total Protein"){myplot1 <- TotalProCheck(exprDat_norm(),bloodStain(),metadata_reord())
    myplot1}
    
    if(input$Assessment == "Calibrator Check"){myplot1 <- CalibratorCheckPlot(calib_norm(),calibIDs())
    myplot1}
    
    if(input$Assessment =="Coefficient Variance Breakdown"){myplot1 <- CVbreakPlot(calib_norm(),calibIDs(),calibPlates())
    myplot1}
    
    if(input$Assessment == "Limit of Detection"){myplot1 <- LoDdetection(RawM)
    myplot1}
    
    if(input$Assessment=="Anova per Protein"){ConfounderTable <- ConfouderCheck(exprDat_norm(),plateID(),metadata_reord(),bloodStain(),sampleAge())
    myplot1 <-confounderPlot2(ConfounderTable)
    myplot1}
  })
  
  output$downloadPlot1 <- downloadHandler(
    filename = function(){
      paste(input$Assessment, '.png', sep='')
    },
    content = function(file) {
      png(file)
      myplot1()
      dev.off
    })
  
  output$ptFinal2 <- renderPlot({
    ptPCAfinal <- PCAglob(MetaRaw(),exprDat_norm(),input$PCAx)
    ptPCAfinal
  })
  
  output$downloadPlot2 <- downloadHandler(
    filename = function(){
      paste(input$PCAx, '.png', sep='')
    },
    content = function(file) {
      ggsave(ptPCAfinal(),filename)
    })
  
  ###output tables  
  output$tblFinal <- renderTable({
    
    if(input$Assessment == "External Validation"){
      CorData_norm1 <- ExtVal1(sandwich_master_Ben,exprDat_norm(),toTest1,toTest2,metadata_reord())
      CorData_norm2 <- ExtVal1(sandwich_master_Historic,exprDat_norm(),toTest1,toTest2,metadata_reord())
      CorData = as.table(as.matrix(rbind(CorData_norm1,CorData_norm2)))
      myTable <- CorData
      myTable}
    
    if(input$Assessment == "KNN Batch Effect Test"){BatchEffectM <- KNNtest(exprDat_norm(),MetaRaw())
    myTable <- BatchEffectM
    myTable}},colnames=FALSE)
  
  output$downloadTable <- downloadHandler(
    filename = function() { 
      paste(input$Assessment, '.csv', sep='')},
    content = function(file) {
      write.csv(myTable(), file)
    }
  )
  
  
}


shinyApp(ui = shinyUI, server = shinyServer)
