
idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu' |
  status['most_general'] == 'HC'
sum(idx)

demo.df <- status[idx,][, c('Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis', 'most_general')]
# demo.df <- status[idx,]
dim(demo.df)

demo.df$most_general <- as.character(demo.df$most_general)

demo.df$most_general[demo.df$most_general == 'greyb'] <- 'probable bacterial'
demo.df$most_general[demo.df$most_general == 'greyu'] <- 'unknown'
demo.df$most_general[demo.df$most_general == 'greyv'] <- 'probable viral'
demo.df$most_general[demo.df$most_general == 'HC'] <- 'healthy control'
demo.df$most_general <- as.factor(demo.df$most_general)

demo.df$most_general
colnames(demo.df)
colnames(demo.df) <- c('Age', 'Sex', 'WBC', 'CRP', 'Presentation', 'Diagnostic Group')

rownames(demo.df) <- seq(1, nrow(demo.df))
# demo.df$`Diagnostic Group`

# View(demo.df)
plotly.table <- demo.df
p <- plot_ly(
  type = 'table',
  columnorder = c(0, 1, 2, 3, 4, 5, 6),
  columnwidth = c(25, 20, 20, 20, 20, 100, 50),
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
api_create(p, filename = "clinical_full_table")



