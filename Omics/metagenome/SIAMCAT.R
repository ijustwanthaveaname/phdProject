rm(list = ls())
setwd("~/phdProject/Omics")
library("SIAMCAT")

data("feat_crc_zeller")
data("meta_crc_zeller")
# feature dataframe
feat.crc.zeller[1:10, 1:5]
# meta data
head(meta.crc.zeller)
str(meta.crc.zeller)
# create label with 2 classes include case
label.crc.zeller <- create.label(
    meta = meta.crc.zeller,
    label = "Group", case = "CRC"
)
# create SIAMCAT object
sc.obj <- siamcat(
    feat = feat.crc.zeller,
    label = label.crc.zeller,
    meta = meta.crc.zeller
)
show(sc.obj)

# feature selection
sc.obj <- filter.features(sc.obj,
    filter.method = 'abundance',
    cutoff = 0.001)

# association analysis
sc.obj <- check.associations(sc.obj)

# cofounder analysis
sc.obj <- check.confounders(
    sc.obj,
    fn.plot = 'confounder_plots.pdf',
    meta.in = NULL,
    feature.type = 'filtered'
)

# data normalization
sc.obj <- normalize.features(
    sc.obj,
    norm.method = "log.unit",
    norm.param = list(
        log.n0 = 1e-06,
        n.p = 2,
        norm.margin = 1
    )
)

# create cross-validation parameters
sc.obj <-  create.data.split(
    sc.obj,
    num.folds = 10,
    num.resample = 2
)

# create train model
sc.obj <- train.model(
    sc.obj,
    method = "lasso"
)
model_type(sc.obj)

# access the models
models <- models(sc.obj)
str(models)
str(models[[1]])
