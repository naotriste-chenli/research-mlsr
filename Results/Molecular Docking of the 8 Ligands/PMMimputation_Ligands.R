library(mice)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

df <- read.csv(args[1])
imp <- mice(df, m = 100, maxit = 20, method = 'pmm')

X11()
densityplot(imp, xlab="Ligand Efficiency")

reference <- complete(imp, 1)
protein_list <- unique(reference$ligand)
final_results <- c()
N <- 1

for (i in protein_list) {

    for (n in 1:100) {
    
        x <- complete(imp, n)
        x_mean <- mean(x$ligand_efficiency)
        x_stdev <- sd(x$ligand_efficiency)
        x_normalized <- (x$ligand_efficiency - x_mean)/x_stdev
        x_max <- max(x_normalized)
        specific_protein <- x %>% filter(grepl(i, x$ligand))
        specific_le <- specific_protein$ligand_efficiency
        normalized_specific <- (specific_le - x_mean)/x_stdev
        normalized_mean <- mean(normalized_specific)
        absolute_difference <- abs(x_max - normalized_mean)/sqrt(2)
        
        final_results[[N]] <- data.frame(ligand = i, distance = absolute_difference)
        
        N <- N + 1
    }

}

final_dataset <- bind_rows(final_results)

anova <- aov(distance ~ ligand, data = final_dataset)
summary(anova)

TukeyHSD(anova)

final_output <- aggregate(cbind(distance) ~ ligand, data = final_dataset, FUN = mean)

ggplot(final_output, aes(x = reorder(ligand, distance), y = distance)) + geom_col(fill = "steelblue") + coord_flip() + theme_minimal() + labs(x = "ID", y = "Distance")

Sys.sleep(4000)

