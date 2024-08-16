##### niche comparisons

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMTools)
library(terra)

##### setup data ----------

# environment/range data
envs <- rast(list.files(path = 'envs', pattern = '.tif$', full.names = T))
envs <- terra::subset(envs, c('bio1', 'bio2', 'bio3', 'bio12', 'bio14'))
plot(envs[[1]])

# occurrence points == imported as SpatVector
g.occs <- vect(read.csv('occs_thinned/green_thinned_thin1.csv'), geom = c('long', 'lat'))
y.occs <- vect(read.csv('occs_thinned/yellow_thinned_thin1.csv'), geom = c('long', 'lat'))

head(g.occs)
head(y.occs)

# background points == imported as SpatVector
bg <- read.csv('bg/bg_10000.csv')
bg$X <- NULL

colnames(bg) = c('long', 'lat')
head(bg)

bg <- vect(bg, geom = c('long', 'lat'))

# define enmtools.species
green <- enmtools.species(range = envs[[1]], presence.points = g.occs, background.points = bg, species.name = 'nigroadspersus_green')
yellow <- enmtools.species(range = envs[[1]], presence.points = y.occs, background.points = bg, species.name = 'nigroadspersus_yellow')


##### Niche identity test ----------
id.test.ecospat <- enmtools.ecospat.id(species.1 = green, species.2 = yellow, env = envs, nreps = 1000, 
                                       R = 100, nback = 10000, bg.source = 'points', verbose = T)

print(id.test.ecospat)

# save results
#saveRDS(id.test.ecospat, 'output_niche_anlyses_rds/id.test.ecospat.rds')


##### background test ----------
bg.test.ecospat <- enmtools.ecospat.bg(species.1 = green, species.2 = yellow, env = envs, nreps = 1000,test.type = 'symmetric', 
                                       R = 100, nback = 10000, bg.source = 'points', verbose = T)

print(bg.test.ecospat)

# save results
#saveRDS(bg.test.ecospat, 'output_niche_anlyses_rds/bg.test.ecospat.rds')

##### plot results----------
