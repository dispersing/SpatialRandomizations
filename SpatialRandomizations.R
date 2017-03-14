##################################################
# This .R file runs Monte Carlo experiments for  #
# to asess correlations between groups of seed   #
# dispersing animals, plants, the difference     #
# between the animals and plants, and            #
# environmental variables (syn. with abiotic).   #
# Each comparison is done in a for loop (it      #
# works and I will write it using apply          #
# functions later) in an array and saved         #
# locally. Datasets are removed ordinally when   #
# they are no longer used later in the file.     #
# The steps are as follows:                      #
# 0. set up                                      #
# 1. prepare abiotic data                        #
# 2. prepare difference data                     #
# 3. simulation of diff ~ abiotics               #
# 4. prepare animal data                         #
# 5. simulation of animals ~ abiotics            #
# 6. prepare plant data                          #
# 7. simulation of plants ~ abiotics             #
# 8. simulation of animals ~ plants              #
##################################################

# 0. set up
	# 0.0. load libraries
		library("sp")
		library("raster")
	# 0.1. set working directory
		wd <- "~/Drive/COOPA_Analysis"
		grd.wd <- paste0(wd, "/GRDs")
		save.wd <- paste0(wd, "/Simulations/Data/animal.plant.Rdata")
		na.wd <- paste0(wd, "/GRDs/NorthAmericaGRDs/na.rast.grd")
	# 0.2. list of files for each set of variables
		abiotic.files <- list.files(path = paste0(grd.wd, "/AbioticGRDs"), pattern = ".grd")
		animal.files <- list.files(path = paste0(grd.wd, "/AnimalGRDs"), pattern = ".grd")
		difference.files <- list.files(path = paste0(grd.wd, "/DifferenceGRDs"), pattern = ".grd")
		plant.files <- list.files(path = paste0(grd.wd, "/PlantGRDs"), pattern = ".grd")
	# 0.3. raster parameters
		final.projection <- "+proj=utm +ellps=WGS84 +datum=WGS84"
		final.extent <- extent(c(-164.4167, -54.41667, 25.06667, 73.56667)) 
		final.rows <- 1077
		final.columns <- 1221
		final.raster <- raster(nrow = final.rows, ncol = final.columns, ext = final.extent, crs = final.projection)
	# 0.4 obtain values for terrestrical north america
		na.rast <- raster(x = na.wd)
		not.nas <- which(!is.na(values(na.rast)))
		rm(na.rast); gc()
	# 0.5. set simulation parameters
		set.seed(6174)
		sims <- 100

# 1. prepare abiotic data
	# 1.0. load abiotic rasters
		aet <- raster(x = paste0(grd.wd, "/AbioticGRDs/aet.grd"))
		dem <- raster(x = paste0(grd.wd, "/AbioticGRDs/dem.grd"))
		ppt <- raster(x = paste0(grd.wd, "/AbioticGRDs/ppt.grd"))
		lat <- raster(x = paste0(grd.wd, "/AbioticGRDs/lat.grd"))
	# 1.1. create some useable objects for the analysis
		abiotic.rasts <- list(aet, dem, ppt, lat)
		abiotic.names <- c("aet", "dem", "ppt", "lat")
		names(abiotic.rasts) <- abiotic.names

# 2. prepare difference data
	# 2.0. load difference rasters
		d.mut <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.mut.grd"))
		d.hoard <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.h.grd"))
		d.frug <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.f.grd"))
		d.hoard.rod <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.h.r.grd"))
		d.hoard.bird <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.h.b.grd"))
		d.frug.mamm <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.f.m.grd"))
		d.frug.bird <- raster(x = paste0(grd.wd, "/DifferenceGRDs/diff.f.b.grd"))
	# 2.1. create some useable objects for the analysis
		diff.rasts <- list(d.mut, d.hoard, d.frug, d.hoard.rod, d.hoard.bird, d.frug.mamm, d.frug.bird)
		diff.names <- c("d.mut", "d.hoard", "d.frug", "d.hoard.rod", "d.hoard.bird", "d.frug.mamm", "d.frug.bird")
		names(diff.rasts) <- diff.names

# 3. simulation of diff ~ abiotics
	# 3.0. preallocate space for randomization
		n.diff <- length(diff.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.diff, n.abiotics))
		dimnames(z.arr) <- list(1:(sims + 1), names(diff.rasts), names(abiotic.rasts))
	# 3.1. loop through comparison and save object
		for (i in 1:n.diff) {
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(diff.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.vals <- values(diff.rasts[[i]])[not.nas]
				for (k in 2:(sims + 1)) {
					z.samp <- sample(z.vals)
					z.rast.samp <- final.raster
					values(z.rast.samp)[not.nas] <- z.samp
					z.arr[k, i, j] <- cor(values(z.rast.samp)[not.nas], values(abiotic.rasts[[j]])[not.nas])			
				} # ends sims, k
			} # ends abiotic maps, j
		} # ends difference maps, i
		diff.abiotic <- z.arr
		save(diff.abiotic, file = paste0(wd, "/Simulations/Data/diff.abiotic.Rdata"))
	# 3.2. remove rasters
		rm(list = diff.names); gc()

# 4. prepare animal data
	# 4.0. load difference rasters
		a.mut <- raster(x = paste0(grd.wd, "/AnimalGRDs/all_muts.grd"))
		a.hoard <- raster(x = paste0(grd.wd, "/AnimalGRDs/all_sh.grd"))
		a.frug <- raster(x = paste0(grd.wd, "/AnimalGRDs/all_frug.grd"))
		a.hoard.rod <- raster(x = paste0(grd.wd, "/AnimalGRDs/rodent_sh.grd"))
		a.hoard.bird <- raster(x = paste0(grd.wd, "/AnimalGRDs/bird_sh.grd"))
		a.frug.mamm <- raster(x = paste0(grd.wd, "/AnimalGRDs/mamm_fruit.grd"))
		a.frug.bird <- raster(x = paste0(grd.wd, "/AnimalGRDs/bird_frug.grd"))
	# 4.1. create some useable objects for the analysis
		animal.rasts <- list(a.mut, a.hoard, a.frug, a.hoard.rod, a.hoard.bird, a.frug.mamm, a.frug.bird)
		animal.names <- c("a.mut", "a.hoard", "a.frug", "a.hoard.rod", "a.hoard.bird", "a.frug.mamm", "a.frug.bird")
		names(animal.rasts) <- animal.names

# 5. simulation of animals ~ abiotics
	# 5.0. preallocate space for randomization
		n.animals <- length(animal.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.animals, n.abiotics))
		dimnames(z.arr) <- list(1:(sims), names(animal.rasts), names(abiotic.rasts))
	# 5.1. loop through comparison and save object
		for (i in 1:n.animals) {
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(animal.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.vals <- values(animal.rasts[[i]])[not.nas]
					z.vector <- vector(mode = "numeric", length = sims)
				for (k in 2:(sims + 1)) {
					z.samp <- sample(z.vals)
					z.rast.samp <- final.raster
					values(z.rast.samp)[not.nas] <- z.samp
					z.vector[(k - 1)] <- cor(values(z.rast.samp)[not.nas], values(abiotic.rasts[[j]])[not.nas])
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.vector
			} # ends abiotic maps, j
		} # ends animal maps, i
		animal.abiotic <- z.arr
		save(animal.abiotic, file = paste0(wd, "/Simulations/Data/animal.abiotic.Rdata"))

# 6. prepare plant data
	# 6.0. load difference rasters
		p.mut <- raster(x = paste0(grd.wd, "/PlantGRDs/total.mut.grd"))
		p.hoard <- raster(x = paste0(grd.wd, "/PlantGRDs/total.hoard.grd"))
		p.frug <- raster(x = paste0(grd.wd, "/PlantGRDs/total.frug.grd"))
		p.hoard.rod <- raster(x = paste0(grd.wd, "/PlantGRDs/total.hoard.rodents.grd"))
		p.hoard.bird <- raster(x = paste0(grd.wd, "/PlantGRDs/total.hoard.birds.grd"))
		p.frug.mamm <- raster(x = paste0(grd.wd, "/PlantGRDs/total.frug.mammal.grd"))
		p.frug.bird <- raster(x = paste0(grd.wd, "/PlantGRDs/total.frug.bird.grd"))
	# 6.1. create some useable objects for the analysis
		plant.rasts <- list(p.mut, p.hoard, p.frug, p.hoard.rod, p.hoard.bird, p.frug.mamm, p.frug.bird)
		names(plant.rasts) <- c("p.mut", "p.hoard", "p.frug", "p.hoard.rod", "p.hoard.bird", "p.frug.mamm", "p.frug.bird")

# 7. simulation of plants ~ abiotics
	# 7.0. preallocate space for randomization
		n.plants <- length(plant.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.plants, n.abiotics))
		dimnames(z.arr) <- list(1:(sims + 1), names(plant.rasts), names(abiotic.rasts))
	# 7.1. loop through comparison and save object
		for (i in 1:n.plants) {
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(plant.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.vals <- values(plant.rasts[[i]])[not.nas]
					z.vector <- vector(mode = "numeric", length = sims)
				for (k in 2:(sims + 1)) {
					z.samp <- sample(z.vals)
					z.rast.samp <- final.raster
					values(z.rast.samp)[not.nas] <- z.samp
					z.vector[(k - 1)] <- cor(values(z.rast.samp)[not.nas], values(abiotic.rasts[[j]])[not.nas])
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.vector
			} # ends abiotic maps, j
		} # ends plant maps, i
		plant.abiotic <- z.arr
		save(plant.abiotic, file = paste0(wd, "/Simulations/Data/plant.abiotic.Rdata"))
	# 7.2. remove rasters
		rm(list = abiotic.names); gc()

# 8. simulation of animals ~ plants
	# 8.0. preallocate space for randomization
		n.animals <- length(animal.rasts)
		n.plants <- length(plant.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.animals, n.plants))
		dimnames(z.arr) <- list(1:(sims + 1), names(animal.rasts), names(plant.rasts))
	# 8.1. loop through comparison and save object
		for (i in 1:n.animals) {
			for (j in 1:n.plants) {
				obs.cor <- cor(values(animal.rasts[[i]])[not.nas], values(plant.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.vals <- values(animal.rasts[[i]])[not.nas]
				for (k in 2:(sims + 1)) {
					z.samp <- sample(z.vals)
					z.rast.samp <- final.raster
					values(z.rast.samp)[not.nas] <- z.samp
					z.arr[k, i, j] <- cor(values(z.rast.samp)[not.nas], values(plant.rasts[[j]])[not.nas])			
				} # ends sims, k
			} # ends plant maps, j
		} # ends animal maps, i
		animal.plant <- z.arr
		save(animal.plant, file = save.wd)