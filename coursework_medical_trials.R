###############################################
# Probability and Statistics Coursework R Code
# B723924
# February-March 2024
###############################################

##############
## Part 1
##############

# Load necessary Libraries
library(tidyverse)
library(dplyr)
library(caTools)

# Set the working directory with the datasets in
setwd("<Path-to-your-data>")

# Clean dataset to remove trailing commas
lines <- readLines("covariates.csv")
lines <- sub(",+$", "", lines)
writeLines(lines, "covariates_cleaning.csv")

# Read the csv files
biomarkers_raw <- read.csv("biomarkers.csv")
covariates_raw <- read.csv("covariates_cleaning.csv")

# Get rid of empty rows
biomarkers_raw <- na.omit(biomarkers_raw)

# Get the patientID into a separate column to allow merging of datasets later
biomarkers_raw$PatientID <- sapply(strsplit(biomarkers_raw$Biomarker, "-"), function(x) x[1])

# Convert to tibble and change data types
biomarkers_raw <- as_tibble(biomarkers_raw)
covariates_raw <- as_tibble(covariates_raw)

biomarkers_raw <- biomarkers_raw %>%
  mutate(PatientID = as.integer(PatientID)) %>%
  arrange(PatientID)

# Merge the two datasets together on PatientID
merged_datasets <- left_join(biomarkers_raw, covariates_raw, by = "PatientID")

# Create a category for high VAS and low VAS to test the hypothesis
merged_datasets <- merged_datasets %>%
  mutate(VAS_at_inclusion_category = ifelse(VAS.at.inclusion >= 5, "High", "Low"))

# As there is only VAS for inclusion remove 12 months and 6 weeks biomarker data
result_tibble <- merged_datasets %>%
  filter(grepl("-0weeks$", Biomarker)) %>%
  select(PatientID, VAS_at_inclusion_category, Biomarker, IL.8, VEGF.A, OPG, TGF.beta.1, IL.6, CXCL9, CXCL1, IL.18, CSF.1) %>%
  arrange(PatientID)

# Subset the data into two groups based on VAS_combined
high_VAS_group <- result_tibble %>% filter(VAS_at_inclusion_category == "High")
low_VAS_group <- result_tibble %>% filter(VAS_at_inclusion_category == "Low")

# Create the tibbles where all the t test results will go
t_test_results <- tibble(Biomarker = character(), t_value = numeric(), p_value_orig = numeric(),p_value_bonferroni = numeric())

# Get the biomarkers from the tibble - columns 4 to 12
biomarkers <- colnames(result_tibble)[4:12]

# For each of the biomarkers, run through the t test and output to tibble
for (biomarker in biomarkers){
  t_test_result <- t.test(high_VAS_group[[biomarker]], low_VAS_group[[biomarker]])
  
  result_row <- tibble(
    Biomarker = biomarker,
    t_value = t_test_result$statistic,
    p_value_orig = t_test_result$p.value,
  )
  t_test_results <- bind_rows(t_test_results, result_row)
}

# d part ii - Bonferroni Corrected results
num_tests <- length(biomarkers)

for (biomarker in biomarkers){
  t_test_result_adj <- t.test(high_VAS_group[[biomarker]], low_VAS_group[[biomarker]])
  
  # Apply Bonferroni correction to the p-value
  adjusted_p_value <- p.adjust(t_test_result_adj$p.value, method = "bonferroni", n=num_tests)

  biomarker_index <- which(t_test_results$Biomarker == biomarker)
  t_test_results$p_value_bonferroni[biomarker_index] <- adjusted_p_value
}

# Display in decimal format and not scientific for reading
t_test_results$p_value_orig <- format(t_test_results$p_value_orig, scientific = FALSE)
t_test_results$p_value_bonferroni <- format(t_test_results$p_value_bonferroni, scientific = FALSE)


############################################################
## Part 2
## Predictions of how well patients with medical condition will recover, 
## 12 Month VAS as response variable and biomarker levels (at inclusion) 
## and covariates as explanatory variables
############################################################

## Take the data set for biomarker levels at inclusion
month0_regression_dataset <- merged_datasets %>%
  filter(grepl("-0weeks$", Biomarker)) %>%
  select(Vas.12months, IL.8, VEGF.A, OPG, TGF.beta.1, IL.6, 
         CXCL9, CXCL1, IL.18, CSF.1, Age, Sex..1.male..2.female., 
         Smoker..1.yes..2.no.)

# Set the seed for reproducibility of the training and test set
set.seed(123)
split_size <- 0.8
split <- sample.split(month0_regression_dataset$Vas.12months, SplitRatio = split_size)
train_data <- subset(month0_regression_dataset, split == TRUE)
test_data <- subset(month0_regression_dataset, split == FALSE)

# Output dimensions of data
dim(train_data)
dim(test_data)

# Run a linear regression model on all the training data
model <- lm(Vas.12months ~ ., data = train_data)

# Get summary data
summary_table <- summary(model)
print(summary_table)

# Plot the residuals both in scatter plot and histogram
res <- resid(model)
plot(fitted(model), res, pch = 16, col = "blue", xlab = "Fitted Values", ylab = "Residuals",)
abline(0,0, col="red")
hist(residuals(model), main = "", xlab = "Residuals")

## Make predictions based on the test data
predictions <- predict(model, newdata = test_data)
predicted_data <- cbind(test_data, Predicted_Vas_12months = predictions)

### Scatter plot comparing predicted vs actual values
par(mar = c(4, 4, 2, 2))
plot(predicted_data$Vas.12months, predicted_data$Predicted_Vas_12months,
     xlab = "Actual 12-month VAS", ylab = "Predicted 12-month VAS",
     main="")
     #main = "Scatter Plot: Actual vs Predicted 12-month VAS")
abline(0, 1, col = "red", lty = 2)