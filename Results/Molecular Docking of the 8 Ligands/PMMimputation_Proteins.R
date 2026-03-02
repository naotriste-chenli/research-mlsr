library(mice)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

df <- read.csv(args[1])
imp <- mice(df, m = 100, maxit = 20, method = 'pmm')
fit <- with(imp, lm(ligand_efficiency ~ complexity))

X11()

pdf("densityplot.pdf", width=8, height=6)

densityplot(imp, xlab="Ligand Efficiency")

dev.off()

reference <- complete(imp, 1)
protein_list <- unique(reference$protein)
final_results <- c()
N <- 1

for (i in protein_list) {

    for (n in 1:100) {
    
        x <- complete(imp, n)
        x_mean <- mean(x$ligand_efficiency)
        x_stdev <- sd(x$ligand_efficiency)
        x_normalized <- (x$ligand_efficiency - x_mean)/x_stdev
        x_max <- max(x_normalized)
        specific_protein <- x %>% filter(grepl(i, x$protein))
        specific_le <- specific_protein$ligand_efficiency
        normalized_specific <- (specific_le - x_mean)/x_stdev
        normalized_mean <- mean(normalized_specific)
        absolute_difference <- abs(x_max - normalized_mean)/sqrt(2)
        
        final_results[[N]] <- data.frame(protein = i, distance = absolute_difference)
        
        N <- N + 1
    }

}

final_dataset <- bind_rows(final_results)
final_output <- aggregate(cbind(distance) ~ protein, data = final_dataset, FUN = mean)

pdf("histogram.pdf", width=8, height=6)

ggplot(final_output, aes(x = reorder(protein, distance), y = distance)) + geom_col(fill = "steelblue") + coord_flip() + theme_minimal() + theme(axis.text.y = element_text(size = 4)) + labs(x = "ID", y = "Distance")

dev.off()

pool(fit)
print(final_output)

Sys.sleep(4000)

