##### Bothriechis nigroadspersus Yellow morph ENM

# set random seed
set.seed(555)

# clean working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMeval)
library(terra)
library(dplyr)

##### part 1 ::: load & prep environment data  ----------
envs <- rast(list.files(path = 'envs', pattern = '.tif$', full.names = T))
envs <- terra::subset(envs, c('bio1', 'bio2', 'bio3', 'bio12', 'bio14'))
plot(envs[[1]])


##### part 2 ::: load occurrence points  ----------

# read raw
y.occs <- read.csv('occs/Yellow.csv')

y.occs$species = 'B.nigroadspersus_yellow'
colnames(y.occs) = c('long', 'lat', 'species')
head(y.occs)

# thin
spThin::thin(loc.data = y.occs, lat.col = 'lat', long.col = 'long', spec.col = 'species', thin.par = 5, reps = 1, locs.thinned.list.return = F, 
             write.files = T, max.files = 1, out.dir = 'occs_thinned', out.base = 'yellow_thinned', write.log.file = F, verbose = T)

# import thinned
y.occs <- read.csv('occs_thinned/yellow_thinned_thin1.csv')
head(y.occs)

plot(envs[[1]])
points(y.occs[, c(2,3)], col = 'yellow')


##### part 3 ::: load background data  ----------
bg <- read.csv('bg/bg_10000.csv')
bg$X <- NULL

colnames(bg) = c('long', 'lat')
head(bg)


##### part 4 ::: data partitioning  ----------
occs.z <- cbind(y.occs[, c(2,3)], terra::extract(envs, y.occs[, c(2,3)]))
bg.z <- cbind(bg, terra::extract(envs, bg))

kfold <- get.randomkfold(occs = y.occs[, c(2,3)], bg = bg, k = 10)
evalplot.envSim.hist(sim.type = 'mess', ref.data = 'occs', occs.z = occs.z, bg.z = bg.z,
                     occs.grp = kfold$occs.grp, bg.grp = kfold$bg.grp)


##### part 5 ::: test models  ----------
eval_yellow <- ENMevaluate(taxon.name = 'B.nigroadspersus_yellow',
                           occs = y.occs[, c(2,3)], 
                           envs = raster::stack(envs),
                           bg = bg,
                           tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'),
                                            rm = seq(0.5, 5, by = 0.5)),
                           partitions = 'randomkfold',
                           partition.settings = list(kfolds = 10),
                           algorithm = 'maxent.jar',
                           doClamp = T)


# select optimal models
eval.yellow_res <- eval.results(eval_yellow)

opt.yellow <- eval.yellow_res %>% filter(cbi.val.avg > 0.85) %>%
  filter(or.10p.avg == min(or.10p.avg))

print(opt.yellow)

# look at the prediction
opt.yellow.pred <- eval.predictions(eval_yellow)[[opt.yellow$tune.args]]
plot(opt.yellow.pred)


##### part 6 ::: binary maps  ----------

# calculate 10% presence threshold
calc_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}

# get threshold
yellow_th <- calc_threshold(sdm = opt.yellow.pred, occs = y.occs[, c(2,3)], type = 'p10', binary = F)

# get binary map
yellow_bin <- ecospat::ecospat.binary.model(Pred = rast(opt.yellow.pred), Threshold = minValue(yellow_th))
plot(yellow_bin)


# save model outputs
#saveRDS(eval_yellow, 'output_models_rds/nigroadspersus_yellow.rds')
#write.csv(eval.yellow_res, 'output_models_metrics/yellow_results.csv')
#writeRaster(opt.yellow.pred, 'output_models_preds/yellow_opt_preds.tif')
#writeRaster(yellow_bin, 'output_models_preds/yellow_bin.tif')
