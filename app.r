#install.packages(deSolve)
library(deSolve)
#install.packages(ggplot2)
library(ggplot2)
#install.packages(dplyr)
library(dplyr)
#install.packages(shiny)
library(shiny)
#library(shinylive)
#library(httpuv)

ebola.app <- shinyApp(
    # This creates the User Interface (UI)
    ui <- pageWithSidebar(
            headerPanel("EBOLA TWO PATCH MODEL"),
            sidebarPanel(
            sliderInput("T", "Time range:",
                              min = 0, max = 100, value = c(0,1),
                        ),
            sliderInput("ATR", "Amplitude Transmission Rate:",
                        min = 0, max = 0.5, value = 0.5,
            ),
            sliderInput("ADR", "Amplitude Dispersal Rate:",
                        min = 0, max = 0.5, value = 0.5,
            ),
            sliderInput("AV_TR_A", "Avg. Transmission Rate (Patch A):",
                        min = 0, max = 60, value = 60,
            ),
            sliderInput("AV_TR_B", "Avg. Transmission Rate (Patch B):",
                        min = 0, max = 5, value = 5,
            ),
            sliderInput("AV_DR_AB", "Avg. Dispersal Rate (AB):",
                        min = 0, max = 0.4, value = 0.4,
            ),
            sliderInput("AV_DR_BA", "Avg. Dispersal Rate (AB):",
                        min = 0, max = 0.04, value = 0.04,
            )

    ),
        mainPanel(
          tabsetPanel(
            
            tabPanel("Model",
            #         withMathJax(

             #          helpText("This is a Two-Patch ODE EBOLA Model")
            #         )
            #),
              tabPanel("", plotOutput("plot1"),
                       helpText("")),
              tabPanel("", plotOutput("plot2"),
              helpText(""))

          ))
        )
    ),

    # This creates the 'behind the scenes' code (Server)
    server <- function(input, output) {

        textsize <- 14

#        mytheme <- theme(legend.position = "top",
 #                        axis.title = element_text(size = textsize + 2),
  #                       axis.text = element_text(size = textsize),
   #                      legend.title = element_text(size = textsize + 2),
    #                     legend.text = element_text(size = textsize))

        two_patch_SIR <- function(time, variables, parameters, phi_BA, phi_AB, beta_A, beta_B) {
  with(as.list(c(variables, parameters)), {
    
    dS_A <- (mu_A * N_A) - ((beta_A(time))*(I_A/N_A)*(S_A)) - ((mu_A + phi_AB(time))*S_A) + (phi_BA(time)*S_B) # Correct.
    
    dI_A <- (beta_A(time)*(I_A/N_A)*S_A) - ((mu_A + gamma_A + phi_AB(time)) * I_A) + (phi_BA(time)*I_B) # Correct.
    
    dR_A <- (gamma_A*I_A)-((mu_A + phi_AB(time))*R_A) + (phi_BA(time)*R_B)  # Correct.
    
#-------------------------------------------------------------------#
    dS_B <-(mu_B * N_B) - ((beta_B(time))*(I_B/N_B)*(S_B)) - ((mu_B + phi_BA(time))*S_B) + (phi_AB(time)*S_A) # Correct.
    
    dI_B <-(beta_B(time)*(I_B/N_B)*S_B) - ((mu_B + gamma_B + phi_BA(time))*I_B) + (phi_AB(time)*I_A) # Correct.
    
    dR_B <-(gamma_B*I_B)-((mu_B + phi_BA(time))*R_B) + (phi_AB(time)*R_A) # Correct.
    
    return(list(c(dS_A, dI_A, dR_A, dS_B, dI_B, dR_B)))
  })
}

      output$plot1 <- renderPlot({

        #####################################################
        # Disease Transmission Rate by (t) with Seasonality #
        #####################################################
        
        beta_A<-function(time){
          return(AV_TR_A*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
        }
        
        beta_B<-function(time){
          
          return(AV_TR_B*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
        }
        ##########################
        # Disease Dispersal Rate #  
        ##########################
        
        phi_AB<-function(time){
          return(AV_DR_AB*(1 + (ADR*cos((2*pi*time/1))))) # Correct.
        }
        
        phi_BA<-function(time){
          
          return(AV_DR_BA*(1 + (ADR*(cos((2*pi*time)/1))))) # Correct.
        }
        
        # Set up initial values for the variables
        initial_values <- c(
          S_A = 95, I_A = 5, R_A = 0,   # Initial values for patch A
          S_B = 995, I_B = 5, R_B = 0   # Initial values for patch B
        )
        
        # Define parameters for the system
        parameters_values <- list(
          N_A = 100,
          N_B = 1000,
          mu_A = 0.06,
          mu_B = 0.02,
          gamma_A = 12,
          gamma_B = 12,
          phi_BA = phi_BA,
          phi_AB = phi_AB,
          beta_A = beta_A,
          beta_B = beta_B
       
        )
        
        #########################################
        #            Dynamical parameters       #
        #---------------------------------------#
        
        time <- seq(input$T[1], input$T[2], by = 1/100)
        ATR<-(input$ATR[1])
        ADR<-(input$ADR[1])
        AV_TR_A<-(input$AV_TR_A[1])
        AV_TR_B<-(input$AV_TR_B[1])
        AV_DR_AB<-(input$AV_DR_AB[1])
        AV_DR_BA<-(input$AV_DR_BA[1])
        
        #########################################
        #            TWO PATCH ODE              #
        #---------------------------------------#
                out <- data.frame(ode(
                  y = initial_values,
                  times = time,
                  func = two_patch_SIR,
                  parms = parameters_values
                ))
        
        out %>% ggplot() +
          geom_line(aes(x = time, y = S_A, colour = "aS"), size = 1, alpha = 0.7) +
          geom_line(aes(x = time, y = I_A, colour = "dR"), size = 1, alpha = 0.7) +
          geom_line(aes(x = time, y = R_A , colour = "bE"), size = 1, alpha = 0.7) +
          scale_x_continuous("Time") +
          scale_colour_discrete("Individuals", labels = c("S", "I", "R")) + ggtitle("PATCH A") + theme(legend.position = "right")

    })
      
      output$plot2 <- renderPlot({
        
        #####################################################
        # Disease Transmission Rate by (t) with Seasonality #
        #####################################################
        seasonality=0.5  #
        ##################
        
        beta_A<-function(time){
          return(AV_TR_A*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
        }
        
        beta_B<-function(time){
          
          return(AV_TR_B*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
        }
        ##########################
        # Disease Dispersal Rate #  
        ##########################
        
        phi_AB<-function(time){
          return(AV_DR_AB*(1 + (ADR*cos((2*pi*time/1))))) # Correct.
        }
        
        phi_BA<-function(time){
          
          return(AV_DR_BA*(1 + (ADR*(cos((2*pi*time)/1))))) # Correct.
        }
        
        # Set up initial values for the variables
        initial_values <- c(
          S_A = 95, I_A = 5, R_A = 0,   # Initial values for patch A
          S_B = 995, I_B = 5, R_B = 0   # Initial values for patch B
        )
        
        # Define parameters for the system
        parameters_values <- list(
          N_A = 100,
          N_B = 1000,
          mu_A = 0.06,
          mu_B = 0.02,
          gamma_A = 12,
          gamma_B = 12,
          phi_BA = phi_BA,
          phi_AB = phi_AB,
          beta_A = beta_A,
          beta_B = beta_B
        )
        
        #########################################
        #            Dynamical parameters       #
        #---------------------------------------#
        time <- seq(input$T[1], input$T[2], by = 1/100)
        ATR<-(input$ATR[1])
        ADR<-(input$ADR[1])
        AV_TR_A<-(input$AV_TR_A[1])
        AV_TR_B<-(input$AV_TR_B[1])
        AV_DR_AB<-(input$AV_DR_AB[1])
        AV_DR_BA<-(input$AV_DR_BA[1])
        
        #########################################
        #            TWO PATCH ODE              #
        #---------------------------------------#
        out <- data.frame(ode(
          y = initial_values,
          times = time,
          func = two_patch_SIR,
          parms = parameters_values
        ))
        
        out %>% ggplot() +
          geom_line(aes(x = time, y = S_B, colour = "aS"), size = 1, alpha = 0.7) +
          geom_line(aes(x = time, y = I_B, colour = "dR"), size = 1, alpha = 0.7) +
          geom_line(aes(x = time, y = R_B , colour = "bE"), size = 1, alpha = 0.7) +
          scale_x_continuous("Time") + ggtitle("PATCH B")+
          scale_colour_discrete("Individuals", labels = c("S", "I", "R")) + theme(legend.position = "right")
        
      })
    }
)

