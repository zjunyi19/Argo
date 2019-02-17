library(shiny)
library(ggplot2)
library(fda)
library(Rcpp)
library(shinyjs)
sourceCpp("main.cpp")

x <- seq(0, 1, length.out = 200) 
epsilon <- rnorm(n = 200, mean = 0, sd = .5) # generate random error
y <- sin(10*x) + exp(x) + epsilon
knots <- seq(0, 1, length.out = 20)
temp <- data.frame(x, y)
M=4
basis <- create.bspline.basis(breaks = knots, norder = M)

ui <- fluidPage(
  titlePanel("B-Spline"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "Choose the log lambda:", min = -10, max = 10, step=0.01, value=3),
      radioButtons("options", "Operation Options: ", 
                   c("Add points" = "add",
                     "Remove points" = "remove",
                     "Toggle points" = "toggle",
                     "Do nothing" = "nothing"),
                   selected = "nothing"),
      actionButton("exclude_toggle", "Toggle points"),
      actionButton("exclude_reset", "Reset")
    ),
    mainPanel(
      plotOutput("plot", height = 300,
                 click = "plot_click",
                 brush = brushOpts(
                   id = "plot_brush"
                 )
      )
    )
  ),
  column(width = 6,
         h4("Brushed points"),
         verbatimTextOutput("brush_info")
  ),
  column(width = 6,
         h4("The coefficients for smoothing splines"),
         verbatimTextOutput("spline")
  )
)

server <- function(input, output) {
  B <- reactiveVal(0)
  sol <- reactiveVal(0)
  dataframe <- reactiveValues()
  dataframe$DT <- temp
  vals <- reactiveValues()
  vals$keeprows <- rep(TRUE, nrow(temp))
  # add row when click add button
  observeEvent({
    input$options == "add"
    input$plot_click
  },{
    add_row <- data.frame(x = input$plot_click$x, y = input$plot_click$y)
    dataframe$DF <- rbind(dataframe$DF, add_row)
    vals$keeprows <- rep(TRUE, nrow(dataframe$DT))
  })
  # remove row when click remove button
  observeEvent(input$rem_point, {
    rem_row <- values$dataframe$DF[-nrow(values$dataframe$DF), ]
    values$dataframe$DF <- rem_row
  })
  
  # Toggle points that are clicked when add button is not clicked
  observeEvent({
    !input$addpoint
    input$plot_click
  }, {
    res <- nearPoints(dataframe$DF, input$plot_click, allRows = TRUE)
    vals$keeprows <- xor(vals$keeprows, res$selected_)
  })
  
  # Toggle points that are brushed when button is clicked
  observeEvent(input$exclude_toggle, {
    res <- brushedPoints(dataframe$DF, input$plot_brush, allRows = TRUE)
    vals$keeprows <- xor(vals$keeprows, res$selected)
  })
  
  # Reset all points
  observeEvent(input$exclude_reset, {
    vals$keeprows <- rep(TRUE, nrow(dataframe$DF))
  })
  
  output$plot <- renderPlot({
    keep <- dataframe$DF[vals$keeprows, ,drop=FALSE]
    exclude <- dataframe$DF[!vals$keeprows, ,drop=FALSE]
    select_x <- brushedPoints(keep, input$plot_brush)$x 
    select_y <- brushedPoints(keep, input$plot_brush)$y
    if (length(select_x) > 1) {
      y_hat <- B() %*% sol()
      predict <- data.frame(select_x, y_hat)
      actual <- data.frame(select_x, select_y)
      ggplot(keep, aes(x,y), color="black") + geom_point() + 
        geom_point(data = exclude, shape = 21, fill = NA, color = "black", alpha = 0.25) +
        geom_line(data=predict, aes(select_x, y_hat), color="red", size = 1)
    } else {
      ggplot(keep, aes(x,y)) + geom_point() + 
        geom_point(data = exclude, shape = 21, fill = NA, color = "black", alpha = 0.25) 
    }
  }) 
  
  # Display brush selected x, y values
  output$brush_info <- renderPrint({
    keep <- dataframe$DF[vals$keeprows, ,drop=FALSE]
    brushedPoints(keep, input$plot_brush)
  })
  
  output$spline <- renderPrint({
    keep <- dataframe$DF[vals$keeprows, ,drop=FALSE]
    exclude <- dataframe$DF[!vals$keeprows, ,drop=FALSE]
    select_x <- brushedPoints(keep, input$plot_brush)$x
    select_y <- brushedPoints(keep, input$plot_brush)$y
    if (length(select_x) <= 1 ) {
      "Please select more data."
    } else {
      B(calcB(select_x, knots, M))
      Omega_B <- fda::bsplinepen(basis)
      A <- t(B()) %*% B() + Omega_B * (10 **input$lambda)
      sol(solve(A) %*% t(B()) %*% select_y)
      sol()
    }
  })
  
  
 

}

shinyApp(ui, server)