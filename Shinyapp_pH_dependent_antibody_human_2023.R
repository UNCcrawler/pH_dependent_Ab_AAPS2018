###########################################################
# Project: Shiny app for simulation of pH-dependent antibody
# Author: Dongfen Yuan
# Date: June 11, 2023
###########################################################
library(shiny)
library(RxODE)
library(tidyverse)

# mPBPK model with a nested endosome compartment
ode = "
 # mPBPK parameters for a 70 kg human
 Vp = 2.6			# L, plasma volume
 Visf = 15.6			# L, total interstitial fluid volume
 Fisf = 0.8		 	# fraction of interstitial fluid available for IgG1 distribution
 V1 = 10.14		 	# L, interstitial fluid volume of lumped tight tissue
 V2 = 5.46			 # L, interstitial fluid volume of lumped leaky tissue
 Vlymph = 5.2	 		# L, lymph volume
 L = 0.121		 	 # L/h, lymph flow
 L1 = 0.0399		 	# L/h, convection from plasma to lumped tight tissue
 L2 = 0.081	 		# L/h, convection from plasma to lumped leaky tissue
 Rf1 = 0.945		 	# vascular reflection coefficient of lumped tight tissue, average data in PMID 25146360 
 Rf2 = 0.697		 	# vascular relection coefficient of lumped leaky tissue, average data in PMID 25146360
 Rflymph = 0.2		 	# lymphatic reflection coefficient
 
 # Endosome parameters
 Ve = 0.0035				# L, endosome volume
 CLup = 0.061667			# L/h, non-specific pinocytosis rate
 CLe = 0.23978			#L/h, manually fitted data for CLup = 0.061667 and CLrec = 0.0189125
 CLrec = 0.693/(8/60)*Ve	# L/h recycling clearance
  
 k1on = 0.0008676			#1/(pM*h), antibody-FcRn binding on rate
 k1off = 583.2				#1/h, antibody-FcRn binding off rate
 FcRntotal = 49.8E6			#pM
 
 # Target parameters
 CLp_t = 0.693/Thalf*Vp  # L/h, total target clearance
 CLp_t1 = CLp_t - CLup	# L/h, endosome-independent target clearance
 ksyn = R0*CLp_t  # pmol/h, target synthesis rate

 
 # Total antibody
 Cp_at = Cp_a + Cp_atc

 
 # Total target
  Cp_tt = Cp_t + Cp_atc
  
 # Initial values
 Cp_a(0) = 0
 Cp_t(0) = R0
 Cp_atc(0) = 0
 Ce_a(0) = 0
 Ce_t(0) = 0
 Ce_atc(0) = 0
 FcRnA(0) = 0
 FcRnATC(0) = 0
 FcRn(0) = FcRntotal
 C1_a(0) = 0
 C1_atc(0) = 0
 C2_a(0) = 0
 C2_atc(0) = 0
 Clymph_a(0) = 0
 Clymph_atc(0) = 0
 

 # Plasma
 
 d/dt(Cp_a) = -kon * Cp_a * Cp_t + koff * Cp_atc - (1-Rf1) * L1 * Cp_a/Vp - (1-Rf2) * L2 * Cp_a/Vp  + L * Clymph_a/Vp - CLup*Cp_a/Vp +CLrec*FcRnA/Vp
 
 d/dt(Cp_t) = ksyn/Vp - CLp_t1 * Cp_t/Vp - CLup*Cp_t/Vp - kon * Cp_a * Cp_t + koff * Cp_atc
 
 d/dt(Cp_atc) = kon * Cp_a * Cp_t - koff * Cp_atc - (1-Rf1) * L1 * Cp_atc/Vp - (1-Rf2) * L2 * Cp_atc/Vp + L * Clymph_atc/Vp -CLup*Cp_atc/Vp +CLrec*FcRnATC/Vp
 
 # Endosome
 
 d/dt(Ce_a) = - keon*Ce_a*Ce_t + keoff*Ce_atc - k1on*FcRn*Ce_a + k1off*FcRnA + CLup*Cp_a/Ve - CLe*Ce_a/Ve
 
 d/dt(Ce_t) = - keon*Ce_a*Ce_t + keoff*Ce_atc - keon*FcRnA*Ce_t + keoff*FcRnATC + CLup*Cp_t/Ve - CLe * Ce_t/Ve
 
 d/dt(Ce_atc) = keon*Ce_a*Ce_t - keoff*Ce_atc - k1on*FcRn*Ce_atc + k1off*FcRnATC + CLup*Cp_atc/Ve - CLe*Ce_atc/Ve
 
 d/dt(FcRnA) = k1on*FcRn* Ce_a - k1off*FcRnA - keon*FcRnA*Ce_t + keoff*FcRnATC - CLrec * FcRnA/Ve
 
 d/dt(FcRnATC) = k1on*FcRn*Ce_atc - k1off*FcRnATC + keon*FcRnA*Ce_t - keoff*FcRnATC - CLrec*FcRnATC/Ve 
 
 d/dt(FcRn) = - k1on*FcRn*Ce_a + k1off*FcRnA - k1on*FcRn*Ce_atc + k1off*FcRnATC + CLrec*(FcRnA + FcRnATC)/Ve
 
 
 # Tight tissue
 
 d/dt(C1_a) = (1-Rf1) * L1 * Cp_a / (V1 * Fisf) - (1-Rflymph) * L1 * C1_a / (V1*Fisf) 
 
 d/dt(C1_atc) = (1-Rf1) * L1 * Cp_atc / (V1 * Fisf) -(1-Rflymph) * L1 * C1_atc / (V1 * Fisf)
 
 # Leaky tissue
 
 d/dt(C2_a) = (1-Rf2) * L2 * Cp_a / (V2 * Fisf) - (1-Rflymph) * L2 * C2_a / (V2 * Fisf)
 
 d/dt(C2_atc) = (1-Rf2) * L2 * Cp_atc / (V2 * Fisf) -(1-Rflymph) * L2 * C2_atc / (V2 * Fisf) 
 
 # Lymph
 
 d/dt(Clymph_a) = (1-Rflymph) * L1 * C1_a / Vlymph + (1-Rflymph) * L2 * C2_a / Vlymph - L * Clymph_a / Vlymph
 
 d/dt(Clymph_atc) = (1-Rflymph) * L1 * C1_atc / Vlymph + (1-Rflymph) * L2 * C2_atc / Vlymph - L * Clymph_atc / Vlymph"


mod = RxODE(model = ode)



ui = fluidPage(

  # App title ----
  titlePanel("Simulation of pH-dependent antibody in human"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Heading
      h4("Antibody"),
      
      # Input: antibody dose ----
      textInput(inputId = "dose",
                   label = "Antibody IV dose (mg/kg, multiple doses separated by comma)",
                   value = "0.25, 0.5, 1, 3, 5"),
      
      # Input: simulation duration ----
      numericInput(inputId = "duration",
                   label = "Simulation duration (days)",
                   value = 30),
      
      br(),
      # Heading
      h4("Soluble target kinetics"),
      
      # Input: soluble target baseline ----
      numericInput(inputId = "R0",
                   label = "Baseline (pM):",
                   value = 0.276),
      
      # Input: soluble target half-life
      numericInput(inputId = "Thalf",
                   label = "Soluble target Half-life (min):",
                   value = 13),
      
      br(),
      # Heading
      h4("Antibody-soluble target binding affinity in plasma"),
      
      # Input: Kon ----
      numericInput(inputId = "Kon",
                   label = "Kon (1/pM*h)",
                   value = 0.002592),
      
      # Input: Koff ----
      numericInput(inputId = "Koff",
                   label = "Koff (1/h):",
                   value = 0.468),
      
      br(),
      # Heading
      h4("Antibody-soluble binding affinity in acidic endosome"),

      # Input: Keon ----
      numericInput(inputId = "Keon",
                   label = "Keon (1/pM*h):",
                   value = 0.002592),
      
      # Input: Keoff ----
      numericInput(inputId = "Keoff",
                   label = "Keoff (1/h):",
                   value = 0.468),
      
      br(),
      # Heading
      h4("Upload observed data"),
      
      # Input: select observed data ----
      fileInput(inputId = "observed",
                label = "Choose .csv file",
                accept = ".csv")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(
        tabPanel(title = "Simulation",
                 # Heading:
                 h4("Antibody plasma PK"),
                 
                 # Output: Antibody PK ----
                 plotOutput(outputId = "PKPlot"),
                 
                 # Heading:
                 h4("Soluble target plasma PK"),
                 
                 # Output: Soluble target PK ----
                 plotOutput(outputId = "PDPlot")
                 
                 ),
        
        tabPanel(title = "About",
                 
                 # Title
                 fixedRow(
                   column(10,
                          h4("A Minimal Physiologically Based Pharmacokinetic Model with a Nested Endosome Compartment for Novel Engineered Antibodies", align = "center"), offset = 0),
                   column(10,
                          h4("Dongfen Yuan, Frederik Rode, Yanguang Cao", align = "center"), offset = 0),
                   column(10,
                          h5("AAPS J. 2018 Mar 14;20(3):48", align = "center"), offset = 0)
                   ),
                 
                 # Section: introduction
                 h4("Introduction"),
                 p("Over the years, several people reached out to me asking questions on this publication. While I am greatful for the interest, I regret that I did not include the model code in the publication. "),
                 p("I took one day and developed a R shiny app. The default setting can be used to reproduce Figure 2b, adalimumab PK in a 70 kg human (digitized data were provided in the data folder. The users can also use it to simulate the PK and PD of pH-dependent antibodies against soluble targets by providing the target baseline/half-life and antibody binding affinities with target in plasma/acidic endosomes."),
                 p("Please refer to the publication for details on model assumptions, methods, and applications."),
                 p("Given that most questions I received are from industry, I chose to upload the source code to Github so that you can download and adapt to your data without risking uploading sensitive data to an external server."),
                 p("If you find errors or have additional questions, please contact Dongfen Yuan via email dongfen.yuan@gmail.com."),
                 p("Dongfen on 11 June 2023"),
                 
                 # Section: model structure
                 h4("Model structure"),
                 
                 img(src = "Model_structure.png", height = 280, width = 700)
                 
                 )
      )
      

      
    )
  )
)




server = function(input, output){
  # RxODE simulation

  sim = reactive({
    theta = c(kon = input$Kon, koff = input$Koff, # antibody-target binding affinity in plasma
              keon = input$Keon, keoff = input$Keoff, # antibody-target binding affinity in acidic endosome
              Thalf = input$Thalf/60, R0 = input$R0)
    
    map_df(strsplit(input$dose, split = ",")%>%unlist%>%as.numeric, 
           function(x){
      ev = eventTable()%>%
        add.dosing(dose = x*70/148000*1000000000/2.6)%>%  # assume mAb 148kDa, mg to pmol/L
        add.sampling(time.units = "h", time = 0:(input$duration*24))
      
      sim = solve(mod, theta, ev)%>%
        select(time, Cp_a, Cp_at, Cp_t, Cp_tt)%>%
        mutate(dose = paste(x, "mg/kg"))
    }) ->simdat
    
    simdat = simdat%>%
      pivot_longer(cols = starts_with("Cp"))%>%
      mutate(time = time/24, # h to day
             value = if_else(name %in% c("Cp_a", "Cp_at"), 
                             value/1000, # pM to nM for antibody
                             value))
    
    return(simdat)
    
  })
  
  obs = reactive({
    if (!is.null(input$observed)){
      
        read.csv(input$observed$datapath, stringsAsFactors = F)%>%
        pivot_longer(starts_with("Cp"))%>%
        mutate(dose = paste(dose, "mg/kg"))
    }
  })
  
  labelname = list(
    Cp_a  = "Free antibody",
    Cp_at = "Total antibody",
    Cp_t = "Free target",
    Cp_tt = "Total target"
  )
  
  # Output: Antibody PK ----
  output$PKPlot = renderPlot({
    
   ggplot()+
      geom_line(data = filter(sim(), name %in% c("Cp_a", "Cp_at")),
                aes(x = time, y = value, color = dose),
                size = 1)+
      facet_wrap(~name, 
                 labeller = function(variable, value){return(labelname[value])})+
      labs(x = "Days", y = "Plasma conc (nM)")+
      scale_y_log10()+
      theme_bw()+
      theme(strip.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))->p
    
    if(!is.null(obs())){
      if(!is.null(filter(obs(), name %in% c("Cp_a", "Cp_at")))){
        p+geom_point(data = filter(obs(), name %in% c("Cp_a", "Cp_at")),
                     aes(x = time/24, y = value, color = dose),
                     size = 3)
      } else{
        p}
    } else {
      p
    }
    

  })
  
  # Output: target PK ----
  output$PDPlot = renderPlot({
    
    ggplot()+
      geom_line(data = filter(sim(), name %in% c("Cp_t", "Cp_tt")),
                aes(x = time, y = value, color = dose),
                size = 1)+
      facet_wrap(~name, 
                 labeller = function(variable, value){return(labelname[value])})+
      labs(x = "Days", y = "Plasma conc (pM)")+
      scale_y_log10()+
      theme_bw()+
      theme(strip.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))-> p
    
    if(!is.null(obs())){
      
      if(!is.null(filter(obs(), name %in% c("Cp_t", "Cp_tt")))){
        p+geom_point(data = filter(obs(), name %in% c("Cp_t", "Cp_tt")),
                     aes(x = time/24, y = value, color = dose),
                     size = 3)
      } else {
        p
      }
    } else {
      p
    }

  })
  
}

shinyApp(ui = ui, server = server)

# img(src = "rstudio.png", height = 140, width = 400)

