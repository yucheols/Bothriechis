##### Bothriechis nigroadspersus Green morph ENM

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
g.occs <- read.csv('occs/Green.csv')

g.occs$species = 'B.nigroadspersus_green'
colnames(g.occs) = c('long', 'lat', 'species')
head(g.occs)

# thin
spThin::thin(loc.data = g.occs, lat.col = 'lat', long.col = 'long', spec.col = 'species', thin.par = 5, reps = 1, locs.thinned.list.return = F, 
             write.files = T, max.files = 1, out.dir = 'occs_thinned', out.base = 'green_thinned', write.log.file = F, verbose = T)

# import thinned
g.occs <- read.csv('occs_thinned/green_thinned_thin1.csv')
head(g.occs)

plot(envs[[1]])
points(g.occs[, c(2,3)], col = 'green')


##### part 3 ::: load background data  ----------
bg <- read.csv('bg/bg_10000.csv')
bg$X <- NULL

colnames(bg) = c('long', 'lat')
head(bg)


##### part 4 ::: data partitioning  ----------
occs.z <- cbind(g.occs[, c(2,3)], terra::extract(envs, g.occs[, c(2,3)]))
bg.z <- cbind(bg, terra::extract(envs, bg))

kfold <- get.randomkfold(occs = g.occs[, c(2,3)], bg = bg, k = 10)
evalplot.envSim.hist(sim.type = 'mess', ref.data = 'occs', occs.z = occs.z, bg.z = bg.z,
                     occs.grp = kfold$occs.grp, bg.grp = kfold$bg.grp)


##### part 5 ::: test models  ----------
eval_green <- ENMevaluate(taxon.name = 'B.nigroadspersus_green',
                          occs = g.occs[, c(2,3)], 
                          envs = raster::stack(envs),
                          bg = bg,
                          tune.args = list(fc = c('L', 'Q', 'H', 'P', 'LQ', 'LP', 'QH', 'QP', 'HP', 'LQH', 'LQP', 'LQHP', 'LQHPT'),
                                           rm = seq(0.5, 5, by = 0.5)),
                          partitions = 'randomkfold',
                          partition.settings = list(kfolds = 10),
                          algorithm = 'maxent.jar',
                          doClamp = T)
