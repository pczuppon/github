library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Burst model"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input values
      numericInput(inputId = "R0", 
                   label = "R0", 
                   value = 13.2)
      ,
      numericInput(inputId = "B", 
                   label = "Burst size (B)", 
                   value = 17.2)
      ,
      numericInput(inputId = "k", 
                   label = "Eclipse rate (k)", 
                   value = 3)
      ,
      numericInput(inputId = "c", 
                   label = "Clearance rate (c)", 
                   value = 10)
      ,
      numericInput(inputId = "T0", 
                   label = "Abundance of target cells (T0)", 
                   value = 4*10^8)
      ,
      # Action button to run
      actionButton(inputId = "clicks", label = "Run")
      ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plot ----
      plotOutput(outputId = "thePlot")
      
    )
  )
)


# SERVER -----------------------------------------
server <- function(input, output) {
  
  # Define generating function minus z to find fixed point
  genfunc <- function(z, parms, epsilon.B, epsilon.bet){
    with(as.list(parms), {
      out <- z - (c/(c + (1 - epsilon.bet)*bet * t) + (1 - epsilon.bet)*bet * t /(c + (1 - epsilon.bet)*bet * t) * exp((1 - epsilon.B)*B * (z-1)))
      return(out)
    })
  }
  
  # Define values of epsilon to be tested
  epsVals <- seq(0, 1, by = 0.001)
  # Initialize output vector
  z.bet <- z.B <- rep(NA, length(epsVals))
  
  # Plot colors
  colB <- "blue"
  colbet <- "orange"

    # If values change, recalculate values
  # Vector of parameters
  parms <- eventReactive(input$clicks, {c(k = input$k, c = input$c, p = input$p, t = input$T0, B = input$B, bet = input$R0 * input$c/(input$T0 * input$B))})
  
  # Plot
  output$thePlot <- renderPlot({
    # Compute the output values by finding roots
    for(i in seq_along(epsVals)){
      tmp <- try(min(1, uniroot(genfunc, c(0, 1-10^(-6)), parms = parms(), epsilon.B = epsVals[i], epsilon.bet = 0, extendInt = "yes")$root), silent = TRUE)
      if(is.numeric(tmp)) z.B[i] <- tmp # NA if no root was found
      
      tmp <- try(min(1, uniroot(genfunc, c(0, 1-10^(-6)), parms = parms(), epsilon.B = 0, epsilon.bet = epsVals[i], extendInt = "yes")$root), silent = TRUE)
      if(is.numeric(tmp)) z.bet[i] <- tmp
    }
    
    # The plot itself
    par(las = 1)
    thelwd <- 2 # line width
    plot(epsVals, 1 - z.B, ylim = c(0,1), type = "l", col = colB, lwd = thelwd,
         xlab = expression(paste("Antiviral drug efficiency ", epsilon)),
         ylab = expression(paste("Establishment probability ", phi)))
    lines(epsVals, 1 - z.bet, ylim = c(0,1), type = "l", col = colbet, lwd = thelwd)
    
    # Add other V values
    Vvals <- c(10, 100)
    for(iV in seq_along(Vvals)){
      lines(epsVals, 1 - z.bet^Vvals[iV], ylim = c(0,1), type = "l", col = colbet, lty = 1 + iV, lwd = thelwd)
      lines(epsVals, 1 - z.B^Vvals[iV], ylim = c(0,1), type = "l", col = colB, lty = 1 + iV, lwd = thelwd)
    }
    
    # Legend
    legend(0, 0, legend = c("reducing burst size B", expression(paste("reducing infectivity ", beta)), "V = 1", "V = 10", "V = 100"), yjust = 0, lty = c(1, 1, 1, 2, 3), col = c(colB, colbet, rep("black", 3)), lwd = thelwd, bty = "n")
    
    })

}

shinyApp(ui = ui, server = server)
