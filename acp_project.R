## Predicting the anticancer potential of peptides with machine learning models ###

# Riku Sundell

### main libraries ----
library(tidyverse)
library(data.table)

### data downloading ----

# download url and destination tempfile
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00589/Anticancer_Peptides.zip"
dl <- tempfile()

# download datasets to tempfile
download.file(url, dl)

# open zip, read in lines of breast and lung sets
breast_set <- fread(text = gsub(",", "\t", readLines(unzip(dl,"ACPs_Breast_cancer.csv"))))
lung_set <- fread(text = gsub(",", "\t", readLines(unzip(dl,"ACPs_Lung_cancer.csv"))))

# combine the two sets into one set by the sequence, rename some parameters, deselect others
combined_set <- full_join(breast_set,
                          lung_set,
                          by = "sequence") %>% 
  rename(breast_cancer = class.x,
         lung_cancer = class.y) %>% 
  select(-ID.x,
         -ID.y)

# remove extra things
rm(breast_set,
   lung_set,
   dl,
   url)

### data structure ----

str(combined_set)
glimpse(combined_set)

### natural amino acids ----

amino_acids <- tibble(amino_acid = c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", "Glutamate", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"), 
                      three_letter = c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"), 
                      one_letter = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                      characteristics = c("Hydrophobic side chain", "Positively charged side chain", "Polar uncharged side chain", "Negatively charged side chain", "Special (contains sulfur)", "Negatively charged side chain", "Polar uncharged side chain", "Special (aminogroup)", "Positively charged side chain", "Hydrophobic side chain", "Hydrophobic side chain", "Positively charged side chain", "Hydrophobic side chain", "Hydrophobic side chain", "Special (cyclic amine)", "Polar uncharged side chain", "Polar uncharged side chain", "Hydrophobic side chain", "Hydrophobic side chain", "Hydrophobic side chain"))

amino_acids

### activities against breast and lung cancer targets ----

combined_set %>%  
  group_by(breast_cancer) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count)) %>%
  ggplot(aes(x = breast_cancer, y = count, fill = breast_cancer)) + geom_bar(stat = "identity", show.legend = FALSE) + coord_flip() +
  labs(title = "Anticancer activity against breast cancer",
       subtitle = "in vitro in cell assay models and/or in silico models",
       x = "",
       y = "Count") +
  theme_classic()

combined_set %>%  
  group_by(lung_cancer) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count)) %>%
  ggplot(aes(x = lung_cancer, y = count, fill = lung_cancer)) + geom_bar(stat = "identity", show.legend = FALSE) + coord_flip() +
  labs(title = "Anticancer activity against lung cancer",
       subtitle = "in vitro in cell assay models and/or in silico models",
       x = "",
       y = "Count") +
  theme_classic()

combined_set %>% 
  group_by(breast_cancer, lung_cancer) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(breast_cancer, lung_cancer, fill = count, label = count)) + 
  geom_tile(show.legend = FALSE) +
  geom_label() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_classic() +
  labs(title = "Anticancer activity against both breast and lung cancer",
       subtitle = "correlation of peptide activity",
       y = "Activity against lung cancer",
       x = "Activity against breast cancer")

# filter the peptides that are very active against both lung and breast cancer
the_seven <- combined_set %>% 
  filter(breast_cancer == "very active" & lung_cancer == "very active")

the_seven

### Calculating new parameters ----

# if you do not have the Peptides package, install it from CRAN
#install.packages("Peptides")

# load the Peptides package
library(Peptides)

aaComp(the_seven$sequence[5])
aaComp(combined_set$sequence[23])

# define pH-range from pH 0 to pH 12 with 0.5 increments
pHs <- seq(from = 0, to = 12, by = 0.5)

#function to run net charge calculation on sequence over pH-range
charge_pH <- function(sequence){
  charge <- charge(sequence, pH = pHs, pKscale = "Bjellqvist")
  tibble(sequence, pHs, charge)
}

#map function over the seven very active ACPs
charge_over_pH_range <- map_df(the_seven$sequence, charge_pH)

#plot the net charge of very active ACPs over pH-range
charge_over_pH_range %>% ggplot(aes(x = pHs, y= charge, color = sequence)) + 
  geom_line() +
  geom_vline(xintercept = 5, color = "gray", lty = 2) + 
  geom_vline(xintercept = 9, color = "gray", lty = 2) +
  labs(title = "Overall charge of 'very active' peptides over a pH-range",
       y = "Net charge",
       x = "pH") +
  theme_classic()

set.seed(749, sample.kind = "Rounding")

inactive <- combined_set %>% filter(breast_cancer == "inactive - exp" & lung_cancer == "inactive - exp") %>% slice_sample(n = 7) %>% select(sequence)

charge_over_pH_range_inactive <- map_df(inactive$sequence, charge_pH)

#plot the net charge of very active ACPs over pH-range
charge_over_pH_range_inactive %>% ggplot(aes(x = pHs, y= charge, color = sequence)) + 
  geom_line() +
  geom_vline(xintercept = 5, color = "gray", lty = 2) + 
  geom_vline(xintercept = 9, color = "gray", lty = 2) +
  labs(title = "Overall charge of 'experimentally inactive' peptides over a pH-range",
       y = "Net charge",
       x = "pH") +
  theme_classic()

#calculate charge at pH 5, 7 and 9
combined_set <- combined_set %>% mutate(charge_pH5 = charge(sequence, pH = 5, pKscale = "Bjellqvist"),
                                        charge_pH7 = charge(sequence, pH = 7, pKscale = "Bjellqvist"),
                                        charge_pH9 = charge(sequence, pH = 9, pKscale = "Bjellqvist"))

# distribution of net charge at pH 5 for breast cancer active...
combined_set %>% ggplot(aes(charge_pH5, fill = breast_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Charge at pH 5") +
  theme_classic()

# ..and lung cancer active peptides
combined_set %>% ggplot(aes(charge_pH5, fill = lung_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Charge at pH 5") +
  theme_classic()

# distribution of net charge at pH 7 with breast cancer active peptides..
combined_set %>% ggplot(aes(charge_pH7, fill = breast_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Charge at pH 7") +
  theme_classic()

# ...and lung cancer active peptides
combined_set %>% ggplot(aes(charge_pH7, fill = lung_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Charge at pH 7") +
  theme_classic()

# distribution of net charge at pH 9 of breast cancer active peptides...
combined_set %>% ggplot(aes(charge_pH9, fill = breast_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Charge at pH 9") +
  theme_classic()

#...and lung cancer active peptides
combined_set %>% ggplot(aes(charge_pH9, fill = lung_cancer)) + 
  geom_density(position = "stack") +
  labs(title = "Net charge of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Charge at pH 9") +
  theme_classic()

#calculating the Boman index
combined_set <- combined_set %>% mutate(boman = boman(sequence))

combined_set %>% ggplot(aes(boman, fill = breast_cancer)) + 
  geom_density(position = "stack") +
  geom_vline(xintercept = 0, color = "gray", lty = 2) +
  theme_classic() +
  labs(title ="The Boman index of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assay",
       x = "Boman index")

combined_set %>% ggplot(aes(boman, fill = lung_cancer)) + 
  geom_density(position = "stack") +
  geom_vline(xintercept = 0, color = "gray", lty = 2) +
  theme_classic() +
  labs(title ="The Boman index of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assay",
       x = "Boman index")

#calculating isoelectric points for the peptides
combined_set <- combined_set %>% mutate(pI = pI(sequence, pKscale = "Bjellqvist"))

combined_set %>% ggplot(aes(pI, fill = breast_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 10, color = "gray", lty = 2) +
  labs(title = "The isoelectric point of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Isoelectric point (pI)") +
  theme_classic()

combined_set %>% ggplot(aes(pI, fill = lung_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 10, color = "gray", lty = 2) +
  labs(title = "The isoelectric point of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Isoelectric point (pI)") +
  theme_classic()

#calculating hydrophobicity of our peptides
combined_set <- combined_set %>% mutate(H = hydrophobicity(sequence, scale = "Eisenberg"))

combined_set %>% ggplot(aes(H, fill = breast_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 0.5, color = "gray", lty =2) +
  labs(title = "Hydrophobicity of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Hydrophobicity") +
  theme_classic()

combined_set %>% ggplot(aes(H, fill = lung_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 0.5, color = "gray", lty =2) +
  labs(title = "Hydrophobicity of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Hydrophobicity") +
  theme_classic()

# calculate the hydrophobic moment of peptides
combined_set <- combined_set %>% mutate(hmoment = hmoment(sequence, angle = 160))

combined_set %>% ggplot(aes(hmoment, fill = breast_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 0.2, color = "gray", lty = 2) +
  labs(title = "Hydrophobic moment of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Hydrophobic moment") +
  theme_classic()

combined_set %>% ggplot(aes(hmoment, fill = lung_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 0.2, color = "gray", lty = 2) +
  labs(title = "Hydrophobic moment of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Hydrophobic moment") +
  theme_classic()

#calculate the instability indexes for our peptides
combined_set <- combined_set %>% mutate(instability = instaIndex(sequence))

combined_set %>% ggplot(aes(instability, fill = breast_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 40, color = "gray", lty = 2) +
  labs(title = "Instability index of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Instability index") +
  theme_classic()

combined_set %>% ggplot(aes(instability, fill = lung_cancer)) +
  geom_density(position = "stack") +
  geom_vline(xintercept = 40, color = "gray", lty = 2) +
  labs(title = "Instability index of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Instability index") +
  theme_classic()

#calculate sequence lenghts
combined_set <- combined_set %>% mutate(peptide_length = lengthpep(sequence))

combined_set %>% ggplot(aes(peptide_length, fill = breast_cancer)) +
  geom_density(position = "stack") + 
  labs(title = "Length of potential anticancer peptides",
       subtitle = "Activity against breast cancer in vitro assays",
       x = "Peptide length") +
  theme_classic()

combined_set %>% ggplot(aes(peptide_length, fill = lung_cancer)) +
  geom_density(position = "stack") + 
  labs(title = "Length of potential anticancer peptides",
       subtitle = "Activity against lung cancer in vitro assays",
       x = "Peptide length") +
  theme_classic()

# count the number of individual amino acids in each sequence
combined_set <- combined_set %>% 
  mutate(A = str_count(sequence, "A"),
         R = str_count(sequence, "R"),
         N = str_count(sequence, "N"),
         D = str_count(sequence, "D"),
         C = str_count(sequence, "C"),
         Q = str_count(sequence, "Q"),
         E = str_count(sequence, "E"),
         G = str_count(sequence, "G"),
         H = str_count(sequence, "H"),
         I = str_count(sequence, "I"),
         L = str_count(sequence, "L"),
         K = str_count(sequence, "K"),
         M = str_count(sequence, "M"),
         F = str_count(sequence, "F"),
         P = str_count(sequence, "P"),
         S = str_count(sequence, "S"),
         T = str_count(sequence, "T"),
         W = str_count(sequence, "W"),
         Y = str_count(sequence, "Y"),
         V = str_count(sequence, "V"))

#count the ratios of residue type vs peptide length
combined_set <- combined_set %>% 
  mutate(ratio_aliphatic = (A+G+I+L+V)/peptide_length,
         ratio_amide = (N+Q)/peptide_length,
         ratio_anion = (D+E)/peptide_length,
         ratio_aromatic = (H+F+W+Y)/peptide_length,
         ratio_cation = (R+K)/peptide_length,
         ratio_cyclic = P/peptide_length,
         ratio_hydroxylic = (S+T)/peptide_length,
         ratio_thiol = C/peptide_length,
         ratio_thioester = M/peptide_length)

# combining the breast and lung cancer activity parameters as one anticancer_activity parameter
combined_set <- combined_set %>% 
  mutate(anticancer_activity = 
           ifelse(lung_cancer == "very active" & breast_cancer == "very active" | 
                    lung_cancer == "very active" & breast_cancer == "mod. active" | 
                    lung_cancer == "mod. active" & breast_cancer == "very active" | 
                    lung_cancer == "mod. active" & breast_cancer == "mod. active" |
                    lung_cancer == "very active" & is.na(breast_cancer) |
                    lung_cancer == "mod. active" & is.na(breast_cancer) |
                    is.na(lung_cancer) & breast_cancer == "very active" |
                    is.na(lung_cancer) & breast_cancer == "mod. active", 
                  "yes", 
                  "no")
  )

### Methods and analysis ----

## Data preprocessing ----

set.seed(789, sample.kind = "Rounding")

#load the caret package
library(caret)

#remove the breast_cancer and lung_cancer parameters
df <- combined_set %>% select(-sequence, -breast_cancer, -lung_cancer)

#create data partition into test and training sets
test_index <- createDataPartition(df$anticancer_activity, times = 1, p = 0.9, list = FALSE)

train_set <- df[test_index,]
test_set <- df[-test_index,]

## Naive model ----

df %>% group_by(anticancer_activity) %>% summarize(count = n())

set.seed(555, sample.kind = "Rounding")
B <- 10000

picks <- sample(df$anticancer_activity, B, replace = TRUE)
prop.table(table(picks))

## Linear model ----

#10 fold cross-validation repeated 10 times
fit_control <- trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 10)
set.seed(6147, sample.kind = "Rounding")

glm_fit <- train(anticancer_activity ~ .,
                 data = train_set,
                 method = "glm",
                 trControl = fit_control)

glm_fit

plot(varImp(glm_fit), top = 20)

glm_pred <- predict(glm_fit, newdata = test_set)

cm_glm <- confusionMatrix(glm_pred, as.factor(test_set$anticancer_activity))

results <- tibble(model = "glm",
                  accuracy = cm_glm$overall["Accuracy"],
                  precision = cm_glm$byClass["Precision"],
                  specificity = cm_glm$byClass["Specificity"],
                  sensitivity = cm_glm$byClass["Sensitivity"],
                  F1_score = cm_glm$byClass["F1"])

library(glmnet)
library(Matrix)

set.seed(42, sample.kind = "Rounding")

glmnet_fit <- train(anticancer_activity ~ .,
                    data = train_set,
                    method = "glmnet",
                    trControl = fit_control)

plot(glmnet_fit)

plot(varImp(glmnet_fit), top = 20)

glmnet_pred <- predict(glmnet_fit, newdata = test_set)

cm_glmnet <- confusionMatrix(glmnet_pred, as.factor(test_set$anticancer_activity))
results <- results %>% add_row(model = "glmnet",
                               accuracy = cm_glmnet$overall["Accuracy"],
                               precision = cm_glmnet$byClass["Precision"],
                               specificity = cm_glmnet$byClass["Specificity"],
                               sensitivity = cm_glmnet$byClass["Sensitivity"],
                               F1_score = cm_glmnet$byClass["F1"])

## KNN ----

#example of optimizing k for knn
control <- trainControl(method = "cv",
                        number = 10,
                        p = 0.9)
#knn training
set.seed(745, sample.kind = "Rounding")

train_knn <- train(anticancer_activity ~., 
                   data = train_set,
                   method = "knn",
                   tuneGrid = data.frame(k = c(1:50)),
                   trControl = control,
                   preProcess = c("center", "scale"),
                   tuneLength = 20)

plot(train_knn)

knn_pred <- predict(train_knn, newdata = test_set)
cm_knn <- confusionMatrix(knn_pred, as.factor(test_set$anticancer_activity))
results <- results %>% add_row(model = "KNN",
                               accuracy = cm_knn$overall["Accuracy"],
                               precision = cm_knn$byClass["Precision"],
                               specificity = cm_knn$byClass["Specificity"],
                               sensitivity = cm_knn$byClass["Sensitivity"],
                               F1_score = cm_knn$byClass["F1"])

## Random Forest ----

library(randomForest)

#10-fold crossvalidation repeated 3 times
control <- trainControl(method="cv", 
                        number = 10,
                        repeats = 3)
#the number of randomly selected variables we will test is 1, 5, 10, 25, 50 and 100
grid <- data.frame(mtry = c(1, 5, 10, 25, 50, 100))

#now we train the model
set.seed(669, sample.kind = "Rounding")

train_rf <-  train(anticancer_activity ~.,
                   data = train_set,
                   method = "rf", 
                   ntree = 150,
                   trControl = control,
                   tuneGrid = grid,
                   nSamp = 5000)

plot(train_rf)

#variable importance for randomForest
plot(varImp(train_rf), top = 20)

rf_pred <- predict(train_rf, newdata = test_set)
cm_rf <- confusionMatrix(rf_pred, as.factor(test_set$anticancer_activity))
results <- results %>% add_row(model = "rf",
                               accuracy = cm_rf$overall["Accuracy"],
                               precision = cm_rf$byClass["Precision"],
                               specificity = cm_rf$byClass["Specificity"],
                               sensitivity = cm_rf$byClass["Sensitivity"],
                               F1_score = cm_rf$byClass["F1"])


## Support Vector Machine ----

#training the model
set.seed(13, sample.kind = "Rounding")

fit_svm <- train(anticancer_activity ~., 
                 data = train_set,
                 method = "svmLinear",
                 trControl = control,
                 preProcess = c("center", "scale"))

svm_pred <- predict(fit_svm, newdata = test_set)

cm_svm1 <- confusionMatrix(svm_pred, as.factor(test_set$anticancer_activity))

results <- results %>% add_row(model = "svmLinear",
                               accuracy = cm_svm1$overall["Accuracy"],
                               precision = cm_svm1$byClass["Precision"],
                               specificity = cm_svm1$byClass["Specificity"],
                               sensitivity = cm_svm1$byClass["Sensitivity"],
                               F1_score = cm_svm1$byClass["F1"])

#training a linear SVM with tuning parameter
set.seed(25, sample.kind = "Rounding")

fit_svm_2 <- train(anticancer_activity ~., 
                   data = train_set,
                   method = "svmLinear",
                   trControl = control,
                   preProcess = c("center", "scale"),
                   tuneGrid = expand.grid(C = seq(0, 2, length = 20)))

plot(fit_svm_2)

svm_pred2 <- predict(fit_svm_2, newdata = test_set)

cm_svm2 <- confusionMatrix(svm_pred2, as.factor(test_set$anticancer_activity))

results <- results %>% add_row(model = "svmLinear with tuning",
                               accuracy = cm_svm2$overall["Accuracy"],
                               precision = cm_svm2$byClass["Precision"],
                               specificity = cm_svm2$byClass["Specificity"],
                               sensitivity = cm_svm2$byClass["Sensitivity"],
                               F1_score = cm_svm2$byClass["F1"])

#training a non-linear SVM with radial kernel and automatic model tuning
set.seed(1997, sample.kind = "Rounding")

fit_svm_3 <- train(anticancer_activity ~., 
                   data = train_set,
                   method = "svmRadial",
                   trControl = control,
                   preProcess = c("center", "scale"),
                   tuneLength = 10)

plot(fit_svm_3)

svm_pred3 <- predict(fit_svm_3, newdata = test_set)

cm_svm3 <- confusionMatrix(svm_pred3, as.factor(test_set$anticancer_activity))

results <- results %>% add_row(model = "svmRadial with tuning",
                               accuracy = cm_svm3$overall["Accuracy"],
                               precision = cm_svm3$byClass["Precision"],
                               specificity = cm_svm3$byClass["Specificity"],
                               sensitivity = cm_svm3$byClass["Sensitivity"],
                               F1_score = cm_svm3$byClass["F1"])

#training a non-linear SVM with polynomial kernel and automatic model tuning
set.seed(11, sample.kind = "Rounding")

fit_svm_4 <- train(anticancer_activity ~., 
                   data = train_set,
                   method = "svmPoly",
                   trControl = control,
                   preProcess = c("center", "scale"),
                   tuneLength = 4)

plot(fit_svm_4)

svm_pred4 <- predict(fit_svm_4, newdata = test_set)

cm_svm4 <- confusionMatrix(svm_pred4, as.factor(test_set$anticancer_activity))

results <- results %>% add_row(model = "svmPoly with tuning",
                               accuracy = cm_svm4$overall["Accuracy"],
                               precision = cm_svm4$byClass["Precision"],
                               specificity = cm_svm4$byClass["Specificity"],
                               sensitivity = cm_svm4$byClass["Sensitivity"],
                               F1_score = cm_svm4$byClass["F1"])

### Results ----
results





