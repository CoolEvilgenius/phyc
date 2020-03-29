library(shiny)
library(markdown)
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Interactive model for genome-based diagnosis of VISA"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    p(strong("Genes")),
    checkboxInput("SA0017",'walR: Response regulator',value = TRUE),
    checkboxInput("SA0018",'walK: Sensor histidine kinase',value = TRUE),
    checkboxInput("SA0500",'rpoB: DNA-directed RNA polymerase β subunit',value = FALSE),
    checkboxInput("SA0573",'sarA: Accessory regulator ',value = FALSE),  
    checkboxInput("SA0614","graR: Response regulator",value = TRUE ), 
    checkboxInput("SA0615","graS: Sensor histidine kinase",value = TRUE  ), 
    checkboxInput("SA0616","vraF: ABC transporter ATP-binding protein" ), 
    checkboxInput("SA0617","vraG: ABC transporter permease"),  
    checkboxInput("SA0723","clpP: Clp protease proteolytic subunit"),  
    checkboxInput("SA1062","stp1: Protein phosphatase 2C"), 
    checkboxInput("SA1557","ccpA: Catabolite control protein"), 
    checkboxInput("SA1659","prsA: Peptidylprolyl isomerase"),  
    checkboxInput("SA1700","vraR: Response regulator",value = TRUE ), 
    checkboxInput("SA1701","vraS: Sensor histidine kinase",value = TRUE ), 
    checkboxInput("SA1702","yvqF: Uncharacterized protein"  ),
    checkboxInput("SA1844","agrA: Accessory gene regulator A",value = TRUE ),   
    p(strong("")),
    p(strong("Threshold for a VISA positive result")),
    p(strong("")),
    sliderInput("etest","Etest MIC",min = 1, max = 8,value = 3, step = 1),
    sliderInput("PAP","PAP-AUC",min = 0.2, max = 2,value = 0.9, step = 0.01),
    p(strong("")),
    p(strong("About this App")),
    includeMarkdown("./description.txt")
    ),
  
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Etest", plotOutput("SSplots"),plotOutput("etestreg")),
      tabPanel("PAP-AUC", plotOutput("PA.SS.plots"), plotOutput("PA.reg"))
    )    
  )
))


