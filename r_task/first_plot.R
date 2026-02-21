# Создаем векторы
genes <- c("BRCA1", "TP53", "EGFR")
expression <- c(12.5, 45.2, 30.1)
condition <- c("Control", "Treatment", "Treatment")

# Объединяем в data frame
exp_data <- data.frame(genes, expression, condition)

# Выводим структуру таблицы
str(exp_data)

# Строим столбчатую диаграмму и сохраняем в PNG
png(filename = "r_task/expression_plot.png")
barplot(exp_data$expression, names.arg = exp_data$genes,
        col = "skyblue", main = "Gene Expression",
        ylab = "Expression Value")
dev.off()
