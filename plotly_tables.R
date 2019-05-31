
demo.df <- status[idx,][, c('most_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
colnames(demo.df) <- c('Diagnosis', 'Site', 'Age', 'Sex', 'WBC', 'CRP', 'Presentation')

rownames(demo.df) <- seq(1, nrow(demo.df))

# View(demo.df)
plotly.table <- demo.df
p <- plot_ly(
  type = 'table',
  columnorder = c(0, 1, 2, 3, 4, 5, 6, 7),
  columnwidth = c(20, 25, 20, 20, 20, 20, 20, 100),
  header = list(
    values = c("<b>Patients</b>", names(plotly.table)),
    align = c('left', rep('center', ncol(plotly.table))),
    line = list(width = 1, color = 'black'),
    fill = list(color = '#444444'),
    font = list(family = "Arial", size = 14, color = "white")
  ),
  cells = list(
    values = rbind(
      rownames(plotly.table), 
      t(as.matrix(unname(plotly.table)))
    ),
    align = c('left', rep('center', ncol(plotly.table))),
    line = list(color = "black", width = 1),
    fill = list(color = c('#9a9e9d')),
    font = list(family = "Arial", size = 12, color = c("black"))
  ))

p
api_create(p, filename = "table_test")



