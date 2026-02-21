# Считываем CSV
data <- read.csv("r_task/sample_data.csv")

# Выводим среднее значение Score
cat("Mean Score:", mean(data$Score), "\n")

# Максимум Score для группы Treatment
max_treatment <- max(data$Score[data$Group == "Treatment"])
cat("Max Score in Treatment:", max_treatment, "\n")

# Строим boxplot по группе
if (!dir.exists("r_task")) {
  dir.create("r_task")
}

png(filename = "r_task/score_boxplot.png")
boxplot(Score ~ Group, data = data,
        main = "Score Distribution by Group",
        col = c("lightgreen","lightcoral"),
        ylab = "Score")
dev.off()
