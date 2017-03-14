# 0. load libraries, directories, and list of files
	# 0.0. load libraries
		library("geoR")
		library("RandomFields")
		library("alr3")
		library("spdep")
		library("fields")
		library("raster")
		library("viridis")
	# 0.1. set directories
		wd <- "~/Drive/COOPA_Analysis"
		grd.wd <- paste0(wd, "/GRDs")
		simdat.wd <- paste0(wd, "/Simulations/Data/GRFData")
	# 0.2. list of files for each set of variables
		abiotic.files <- list.files(path = paste0(grd.wd, "/AbioticGRDs"), pattern = ".grd")
		animal.files <- list.files(path = paste0(grd.wd, "/AnimalGRDs"), pattern = ".grd")
		difference.files <- list.files(path = paste0(grd.wd, "/DifferenceGRDs"), pattern = ".grd")
		plant.files <- list.files(path = paste0(grd.wd, "/PlantGRDs"), pattern = ".grd")
	# 0.3. raster parameters
		final.projection <- "+proj=utm +ellps=WGS84 +datum=WGS84"
		final.extent <- extent(c(-164.4167, -54.41667, 25.06667, 73.56667)) 
		final.columns <- 100 #1221 # longitude, x
		final.rows <- 100 #1077 # latitude, y
		final.raster <- raster(nrow = final.rows, ncol = final.columns, ext = final.extent, crs = final.projection)
	# 0.4. set simulation parameters
		set.seed(6174)
		sims <- 1000

# 1. prepare abiotic data
	# 1.0. load abiotic rasters
		aet <- raster(x = paste0(grd.wd, "/AbioticGRDs/aet.grd"))
		dem <- raster(x = paste0(grd.wd, "/AbioticGRDs/dem.grd"))
		ppt <- raster(x = paste0(grd.wd, "/AbioticGRDs/ppt.grd"))
		lat <- raster(x = paste0(grd.wd, "/AbioticGRDs/lat.grd"))
	# 1.1. create some useable objects for the analysis
		abiotic.rasts.fullres <- list(aet, dem, ppt, lat)
		abiotic.names <- c("aet", "dem", "ppt", "lat")
		names(abiotic.rasts.fullres) <- abiotic.names
	# 1.2. resize rasters
		abiotic.rasts <- lapply(abiotic.rasts.fullres, function(x) resample(x, final.raster, method = "bilinear"))

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
		diff.rasts.fullres <- list(d.mut, d.hoard, d.frug, d.hoard.rod, d.hoard.bird, d.frug.mamm, d.frug.bird)
		diff.names <- c("d.mut", "d.hoard", "d.frug", "d.hoard.rod", "d.hoard.bird", "d.frug.mamm", "d.frug.bird")
		names(diff.rasts.fullres) <- diff.names
	# 2.2. resize rasters
		diff.rasts <- lapply(diff.rasts.fullres, function(x) resample(x, final.raster, method = "bilinear"))

# 3. simulation of diff ~ abiotics
	# 3.0. preallocate space for randomization
		n.diff <- length(diff.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.diff, n.abiotics))
		dimnames(z.arr) <- list(1:(sims + 1), names(diff.rasts), names(abiotic.rasts))
	# 3.1. loop through comparison and save object
		for (i in 1:n.diff) {
				not.nas <- which(!is.na(values(diff.rasts[[i]])) == T)
				z.vals <- values(diff.rasts[[i]])[not.nas]
				dat <- as.data.frame(cbind(coordinates(diff.rasts[[i]]), values(diff.rasts[[i]])))
				colnames(dat) <- c("x", "y", "z")
				gDat <- as.geodata(dat)
				vObs <- variog(gDat)
				VF <- variofit(vObs, cov.model = "gaussian", fix.nugget = F, messages = F, weights = "cressie", max.dist = 80)
				modObj <- geoR2RF(cov.model = VF$cov.model, cov.pars = VF$cov.pars, nugget = VF$nugget, kappa = VF$kappa)
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(diff.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.arr.vec <- vector(mode = "numeric", length = sims)
				for (k in 1:sims) {
					GRF <- RFsimulate(model = modObj, x = 1:100, y = 1:100, grid = TRUE)
					z.arr.vec[k] <- cor(GRF[[1]][not.nas], values(abiotic.rasts[[j]])[not.nas], method = "spearman")
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.arr.vec
			} # ends abiotic maps, j
		} # ends difference maps, i
		diff.abiotic <- z.arr
		save(diff.abiotic, file = paste0(simdat.wd, "/diff.abiotic.Rdata"))
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
		animal.rasts.fullres <- list(a.mut, a.hoard, a.frug, a.hoard.rod, a.hoard.bird, a.frug.mamm, a.frug.bird)
		animal.names <- c("a.mut", "a.hoard", "a.frug", "a.hoard.rod", "a.hoard.bird", "a.frug.mamm", "a.frug.bird")
		names(animal.rasts.fullres) <- animal.names
	# 4.2. resize rasters
		animal.rasts <- lapply(animal.rasts.fullres, function(x) resample(x, final.raster, method = "bilinear"))

# 5. simulation of animals ~ abiotics
	# 5.0. preallocate space for randomization
		n.animals <- length(animal.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.animals, n.abiotics))
		dimnames(z.arr) <- list(1:(sims + 1), names(animal.rasts), names(abiotic.rasts))
	# 5.1. loop through comparison and save object
		for (i in 1:n.animals) {
				not.nas <- which(!is.na(values(animal.rasts[[i]])) == T)
				z.vals <- values(animal.rasts[[i]])[not.nas]
				dat <- as.data.frame(cbind(coordinates(animal.rasts[[i]]), values(animal.rasts[[i]])))
				colnames(dat) <- c("x", "y", "z")
				gDat <- as.geodata(dat)
				vObs <- variog(gDat)
				VF <- variofit(vObs, cov.model = "gaussian", fix.nugget = F, messages = F, weights = "cressie", max.dist = 80)
				modObj <- geoR2RF(cov.model = VF$cov.model, cov.pars = VF$cov.pars, nugget = VF$nugget, kappa = VF$kappa)
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(animal.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.arr.vec <- vector(mode = "numeric", length = sims)
				for (k in 1:sims) {
					GRF <- RFsimulate(model = modObj, x = 1:100, y = 1:100, grid = TRUE)
					z.arr.vec[k] <- cor(GRF[[1]][not.nas], values(abiotic.rasts[[j]])[not.nas], method = "spearman")
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.arr.vec
			} # ends abiotic maps, j
		} # ends animal maps, i
		animal.abiotic <- z.arr
		save(animal.abiotic, file = paste0(simdat.wd, "/animal.abiotic.Rdata"))


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
		plant.rasts.fullres <- list(p.mut, p.hoard, p.frug, p.hoard.rod, p.hoard.bird, p.frug.mamm, p.frug.bird)
		names(plant.rasts.fullres) <- c("p.mut", "p.hoard", "p.frug", "p.hoard.rod", "p.hoard.bird", "p.frug.mamm", "p.frug.bird")
	# 6.2. resize rasters
		plant.rasts <- lapply(plant.rasts.fullres, function(x) resample(x, final.raster, method = "bilinear"))

# 7. simulation of plants ~ abiotics
	# 7.0. preallocate space for randomization
		n.plants <- length(plant.rasts)
		n.abiotics <- length(abiotic.rasts)
		z.arr <- array(data = NA, dim = c((sims + 1), n.plants, n.abiotics))
		dimnames(z.arr) <- list(1:(sims + 1), names(plant.rasts), names(abiotic.rasts))
	# 7.1. loop through comparison and save object
		for (i in 1:n.plants) {
				not.nas <- which(!is.na(values(plant.rasts[[i]])) == T)
				z.vals <- values(plant.rasts[[i]])[not.nas]
				dat <- as.data.frame(cbind(coordinates(plant.rasts[[i]]), values(plant.rasts[[i]])))
				colnames(dat) <- c("x", "y", "z")
				gDat <- as.geodata(dat)
				vObs <- variog(gDat)
				VF <- variofit(vObs, cov.model = "gaussian", fix.nugget = F, messages = F, weights = "cressie", max.dist = 80)
				modObj <- geoR2RF(cov.model = VF$cov.model, cov.pars = VF$cov.pars, nugget = VF$nugget, kappa = VF$kappa)
			for (j in 1:n.abiotics) {
				obs.cor <- cor(values(plant.rasts[[i]])[not.nas], values(abiotic.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.arr.vec <- vector(mode = "numeric", length = sims)
				for (k in 1:sims) {
					GRF <- RFsimulate(model = modObj, x = 1:100, y = 1:100, grid = TRUE)
					z.arr.vec[k] <- cor(GRF[[1]][not.nas], values(abiotic.rasts[[j]])[not.nas], method = "spearman")
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.arr.vec
			} # ends abiotic maps, j
		} # ends plant maps, i
			plant.abiotic <- z.arr
			save(plant.abiotic, file = paste0(simdat.wd, "/plant.abiotic.Rdata"))
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
				not.nas <- which(!is.na(values(animal.rasts[[i]])) == T)
				z.vals <- values(plant.rasts[[i]])[not.nas]
				dat <- as.data.frame(cbind(coordinates(animal.rasts[[i]]), values(animal.rasts[[i]])))
				colnames(dat) <- c("x", "y", "z")
				gDat <- as.geodata(dat)
				vObs <- variog(gDat)
				VF <- variofit(vObs, cov.model = "gaussian", fix.nugget = F, messages = F, weights = "cressie", max.dist = 80)
				modObj <- geoR2RF(cov.model = VF$cov.model, cov.pars = VF$cov.pars, nugget = VF$nugget, kappa = VF$kappa)
			for (j in 1:n.plants) {
				obs.cor <- cor(values(animal.rasts[[i]])[not.nas], values(plant.rasts[[j]])[not.nas])
				z.arr[1, i, j] <- obs.cor
				z.arr.vec <- vector(mode = "numeric", length = sims)
				for (k in 1:sims) {
					GRF <- RFsimulate(model = modObj, x = 1:100, y = 1:100, grid = TRUE)
					z.arr.vec[k] <- cor(GRF[[1]][not.nas], values(plant.rasts[[j]])[not.nas], method = "spearman")
				} # ends sims, k
				z.arr[2:(sims + 1), i, j] <- z.arr.vec
			} # ends plant maps, j
		} # ends animal maps, i
		animal.plant <- z.arr
		save(animal.plant, file = paste0(simdat.wd, "/animal.plant.Rdata"))