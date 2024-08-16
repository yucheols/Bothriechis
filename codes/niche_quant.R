##### niche comparisons

# clear working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ENMTools)
library(terra)
library(ggplot2)
library(ggpubr)
library(dplyr)

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


##### plot results ----------

# id test
print(id.test.ecospat$test.results$sim$D)

sim.D <- as.data.frame(id.test.ecospat$test.results$sim$D)
sim.D$test_type = 'identity'
colnames(sim.D) = c('D', 'test_type')
head(sim.D)

idplot <- sim.D %>%
  ggplot(aes(x = D, fill = test_type, color = test_type)) +
  geom_histogram(bins = 10, alpha = 0.4, linewidth = 1.0) +
  geom_vline(xintercept = id.test.ecospat$test.results$obs$D, color = 'black', linewidth = 1.0, linetype = 'longdash') +
  scale_fill_manual(values = 'cornflowerblue') +
  scale_color_manual(values = 'cornflowerblue') +
  xlab("Schoener's D") + ylab('Count') +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = margin(l = 15, r = 15)) 


# bg test
print(bg.test.ecospat$test.results$sim$D)

bg.D <- as.data.frame(bg.test.ecospat$test.results$sim$D)
bg.D$test_type = 'background'
colnames(bg.D) = c('D', 'test_type')
head(bg.D)

bgplot <- bg.D %>%
  ggplot(aes(x = D, fill = test_type, color = test_type)) +
  geom_histogram(bins = 10, alpha = 0.4, linewidth = 1.0) +
  geom_vline(xintercept = bg.test.ecospat$test.results$obs$D, color = 'black', linewidth = 1.0, linetype = 'longdash') +
  scale_fill_manual(values = 'salmon') +
  scale_color_manual(values = 'salmon') +
  xlab("Schoener's D") + ylab('Count') +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = margin(l = 15, r = 15)) 


# merge plots
ggarrange(idplot, bgplot,
          font.label = list(size = 18),
          labels = c('A)', 'B)'),
          nrow = 1, ncol = 2)

ggsave('plots/niche_analyses_Espace.png', width = 20, height = 10, dpi = 800, units = 'cm')
