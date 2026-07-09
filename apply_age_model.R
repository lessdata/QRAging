# Age prediction and uncertainty quantification using quantile regression (QR)
# Models were trained using UK Biobank Olink proteomics data
# Reference (TBA)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(splines2)
# library(ggplot2)
# library(survival)
options(scipen = 999)

dir = paste0(".")

## 1. Phenotype data
## 1.1 Chronological age at sample collection and other covarites such as sex
df.age = fread(file = paste0(dir, "/example/example.age.tsv"),
               header = T, sep = "\t", data.table = F, stringsAsFactors = F)

## 1.2 Protein expression matrix
## IID column: unique sample identifier
## remaining columns: protein expression values, with one column per protein and one row per sample
df.prot = fread(file = paste0(dir, "/example/example.proteomics.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
mt.prot = df.prot %>% select(-IID) %>% scale() %>% as.data.frame()
mt.prot[is.na(mt.prot)] = 0
mt.prot = mt.prot %>% mutate(Intercept = 1)

## 2. Organ age models
## Select the identifier column from the model weight table to 
## match the protein names in expression matrix.
## For example, 
## if proteins are named using the Olink panel name and gene symbol,
## create an identifier column in the model weight table: 
## df.beta = mutate(Predictor = paste0(Protein.Panel, "_", Gene.Symbol))
## By default, 
## the assay target ID is assumed to be the common identifier
## used to map the proteins between the expression matrix and the model weights
df.beta = fread(file = paste0(dir, "/model/qr.beta.proteomic.age.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
df.beta = df.beta %>% mutate(Predictor = Assay.Target)

## 3. Organ age predictions
df.pred = data.frame()
organ.pool = df.beta %>% distinct(Organ) %>% pull(Organ)
for (organ in organ.pool) {
  #### 3.1 Applying the QR model to the test set
  beta = df.beta %>% 
    filter(Organ == organ,
           Predictor %in% colnames(mt.prot)) %>%
    select(Predictor, starts_with("QR"))
  pred = as.data.frame(as.matrix(mt.prot %>% select(all_of(beta$Predictor))) %*% as.matrix(beta %>% select(-Predictor))) %>%
    mutate(df.prot %>% select(IID), .before = 1)
  
  #### 3.2 Post hoc sorting
  df.qntl = data.frame(Name = pred %>% select(starts_with("QR")) %>% colnames(), 
                       stringsAsFactors = F) %>%
    mutate(Idx = 1:nrow(.))
  pred = pred %>% 
    pivot_longer(cols = starts_with("QR"), names_to = "Name", values_to = "Value") %>%
    group_by(IID) %>%
    mutate(Idx = row_number()) %>%
    arrange(Value, Idx, .by_group = T) %>%
    mutate(Idx = row_number()) %>%
    ungroup() %>%
    select(-Name) %>%
    left_join(y = df.qntl, by = join_by(Idx)) %>%
    pivot_wider(id_cols = IID, 
                names_from = Name, values_from = Value)
  pred = pred %>% 
    left_join(y = df.age, by = join_by(IID)) %>%
    relocate(Sex, Age, .after = IID)
  
  #### 3.3 Age-calibrated age gap and standardized age gap z-score
  #### age prediction => Age.QR0.5
  #### age-calibrated age gap => Gap.QR0.5
  #### standardized age-calibrated age gap z-score => Gap.Std.QR0.5
  fit = lm(data = pred, formula = QR0.5 ~ Age)
  adj = unname(fit$fitted.values)
  pred = pred %>% 
    mutate(across(.cols = starts_with("QR"),
                  .fns = ~ .x - adj,
                  .names = "Gap.{.col}"),
           across(.cols = starts_with("QR"),
                  .fns = ~ .x - adj,
                  .names = "Gap.Std.{.col}")) %>%
    rename_with(.cols = starts_with("QR"), 
                .fn = ~ paste0("Age.", .x))
  m.pred = pred %>% pull(Gap.Std.QR0.5) %>% mean()
  sd.pred = pred %>% pull(Gap.Std.QR0.5) %>% sd()
  pred = pred %>% 
    mutate(across(.cols = starts_with("Gap.Std"),
                  .fns = ~ (.x - m.pred) / sd.pred))
  
  #### 3.4 Tail probability P(Gap > 0)
  #### tail probability => Prob
  pred = pred %>% mutate(Prob = 0)
  th.prob = 0
  n.sample = 10000
  set.seed(seed = 256)
  prob.int = runif(n.sample, 0, 1)
  prob.cum = c(0.025, (1:9) / 10, 0.975)
  set.seed(seed = 256)
  for (idx in 1:nrow(pred)) {
    tmp = pred %>% 
      slice(idx) %>% 
      select(starts_with("Gap.QR")) %>%
      pivot_longer(cols = starts_with("Gap.QR"), 
                   names_to = "Quantile", 
                   values_to = "QR") %>% 
      mutate(Prob = prob.cum)
    fit.spline = lm(data = tmp, formula = QR ~ bSpline(Prob, Boundary.knots = c(0, 1)))
    pred.spline = predict(object = fit.spline, newdata = data.frame(Prob = prob.int))
    pred$Prob[idx] = sum(pred.spline >= th.prob) / n.sample
  }
  
  df.pred = rbind(df.pred,
                  pred %>% mutate(Organ = organ, .before = 1))
  print(paste0("==== ", organ, " ==== ", Sys.time(), " ===="))
}

#### 3.5 Age gap interval length
df.pred = df.pred %>% mutate(Interval.Length = Gap.QR0.975 - Gap.QR0.025)

#### 3.6 Extreme ager status (aging outlier) by age gap
th.z = 1.5
df.pred = df.pred %>% 
  mutate(Outlier.Gap = case_when(Gap.Std.QR0.5 > th.z ~ "Extreme age",
                                 Gap.Std.QR0.5 < -th.z ~ "Extreme youth"
                                 .default = "Normal"))

## 4. Use age gap, tail probability, and interval length as aging measures
## age gap => Gap.Std.QR0.5
## tail probability => Prob
## interval length => Interval.Length
ofile = paste0(dir, "/example/example.prediction.tsv")
fwrite(x = df.pred,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)

