library(rmarkdown)

# Render flexdashboard to HTML
render(
  input  = "Di-GRAPH_report.Rmd",
  output_format = "flexdashboard::flex_dashboard",
  output_file   = "Di-GRAPH_report.html"
)