# Deconvolution of GTEx Coronary Artery Supplementary Materials

This repository contains supplementary information for the deconvolution of GTEx coronary artery tissues. Additional figures are listed separately in this repository, with the following documentation showing how to access the data used in the main work. This requires the recount3 package from Bioconductor and the three text files in this repository srp.txt, sra.txt, and celltypes.txt.

## Download the GTEx data

We used the recount3 package to download the GTEx coronary artery data. These data are
available as a part of the overall blood vessel set found in recount3.

```{r download gtex data}
library(recount3)
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "BLOOD_VESSEL" & project_type == "data_sources"
)
rse_coronary <- create_rse(proj_info)
rse_coronary <- rse_coronary[, rse_coronary$gtex.smtsd == "Artery - Coronary"]
assay(rse_coronary, "counts") <- transform_counts(rse_coronary)
saveRDS(rse_coronary, file = "data/rse_coronary_recount3_full.Rds")
```

There are several subjects with duplicated samples in the coronary artery data so we
subset these data to only include one sample per subject. Features with counts less than
1,000 in all subjects are removed as well. Numerous sample qualities are added to the column data
to account for adjacent histology and demographic information in the GTEx samples.

```{r filter low counts and duplicate subjects}
assays(rse_coronary)[1] <- NULL
lowCounts <- apply(assay(rse, "counts"), 1, function(x) {
  all(x < 1000)
})

rse_coronary <- rse_coronary[-which(lowCounts), ]

gtexNames <- str_split(cibAbs$InputSample, "-") %>%
  sapply(FUN = function(x) paste(x[1], x[2], x[3], sep = "-"))
cibAbs$InputSample <- gtexNames
rse_coronary$histologyID <- gtexNames

duplicatedSamples <- which(duplicated(rse_coronary$gtex.subjid))

## duplicated samples results 88, 120, 212, 213, 214
## 49 GTEX-P4PP-3026 (same subject for 88 and 120, both vein in histology)
## 22 GTEX-QV44-0726 (same subject for 212 which has NA SRR)
## 96 GTEX-QEG5-1026 (same subject for 213 which has NA SRR)
## 32 GTEX-QLQ7-0326 (same subject for 214 which has NA SRR)

rse_coronary <- rse_coronary[, -duplicatedSamples]

cac <- read_excel("Coronary_artery_case_images_8.27.20.xlsx")

cacSub <- subset(cac, `Specimen ID` %in% rse_coronary$histologyID)
## one duplicated Specimen ID
cacSub <- cacSub[-148, ]

map <- match(rse_coronary$histologyID, cacSub$`Specimen ID`)

cacSub <- cacSub[map, ]

## confirm that these tables are aligned, following line should evalute as TRUE

identical(rse_coronary$histologyID, cacSub$`Specimen ID`)

rse_coronary$virmani <- cacSub$`Plaque Type`

rse_coronary$plaque <- factor(rse_coronary$virmani)

rse_coronary$plaque[!(rse_coronary$virmani == 0)] <- 1

rse_coronary$virmani <- factor(rse_coronary$virmani)

rse_coronary$plaque <- factor(rse_coronary$plaque)

rse_coronary$year <- parse_date_time(rse_coronary$gtex.smgebtchd, "mdy") %>%
  year() %>%
  factor()

rse_coronary$age <- factor(rse_coronary$gtex.age, labels = 0:5)

rse_coronary$sex <- factor(rse_coronary$gtex.sex, labels = c("Male", "Female"))

veinSamples <- c(
  "GTEX-15DDE-1226",
  "GTEX-P4PP-3026",
  "GTEX-132AR-1326",
  "GTEX-15TU5-1226",
  "GTEX-RTLS-0726",
  "GTEX-12ZZX-0626",
  "GTEX-1CB4E-0726",
  "GTEX-15CHR-0526",
  "GTEX-ZYT6-0826",
  "GTEX-12WSC-0926",
  "GTEX-14BMV-0926",
  "GTEX-17GQL-0426",
  "GTEX-ZDYS-0126"
)

veinInd <- which(rse_coronary$histologyID %in% veinSamples)

rse_coronary <- rse_coronary[, -veinInd]

saveRDS(rse_coronary,
  file = "data/rse_coronary_recount3_for_cibersort.Rds"
)
```

## Download the SRA data

```{r}
library(recount3)

# load in the samples of interest, their respective celltypes,
# and the project IDs that they are found in.

projects <- read.table("data/srp.txt")
runs <- read.table("data/srr.txt")
celltypes <- read.table("data/celltypes.txt")

# identify the projects available in recount3.
human_projects <- available_projects()

projects <- projects[projects %in% human_projects$project]

project <- subset(
  human_projects,
  project == projects[1] & project_type == "data_sources"
)

rse <- create_rse(project_info = project)

projects <- projects[-1]

for (i in 1:length(projects)) {
  project <- subset(human_projects, project == projects[i] & project_type == "data_sources")
  tmp_rse <- create_rse(project_info = project)
  rse <- cbind(rse, tmp_rse)
}

rse <- rse[, colnames(rse) %in% runs]

map <- match(runs, colnames(rse))

rse <- rse[, map]

rse$celltype <- celltypes

assay(rse, "counts") <- transform_counts(rse)

saveRDS(rse, file = "data/rse_reference_full.Rds")

duplicated_srs <- unique(rse$sra.sample_acc.x[duplicated(rse$sra.sample_acc.x)])

dup_srs_ind <- sapply(duplicated_srs, function(x) which(rse$sra.sample_acc.x == x))

counts_mat <- assay(rse, "counts")

rep_sums <- lapply(dup_srs_ind, function(x) rowSums(counts_mat[, x]))

dup_remove <- sapply(1:length(dup_srs_ind), function(x) dup_srs_ind[[x]][-1]) %>% unlist()

dup_replace <- sapply(1:length(dup_srs_ind), function(x) dup_srs_ind[[x]][1])

for (i in 1:length(rep_sums)) {
  counts_mat[, dup_replace[i]] <- rep_sums[[i]]
}

counts_mat <- counts_mat[, -dup_remove]

rse <- rse[, -dup_remove]

assay(rse, "counts") <- counts_mat

# raw counts not necessary for analyses, remove this assay to conserve memory

assays(rse)[1] <- NULL

lowCounts <- apply(assay(rse), 1, function(x) {
  all(x < 1000)
})

rse <- rse[-which(lowCounts), ]

saveRDS(rse, file = "data/rse_reference_no_reps.Rds")

rse <- rse[, -which(rse$celltype == "Erythrocyte")]

saveRDS(rse, file = "data/rse_reference_noRBC.Rds")
```
