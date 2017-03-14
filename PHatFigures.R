##################################################
# This .R file imports files from the            #
# Monte Carlo simulations run in                 #
# SpatialRansomizations.R.  It loops over the    #
# array objects, plots multi-panel figures, and  #
# saves them as .pdfs to a "/Figures" directory. #
# The steps are as follows:                      #
# 0. set up                                      #
# 1. plot difference ~ abiotic                   #
# 2. plot animal ~ abiotic                       #
# 3. plot plant ~ abiotic                        #
# 4. plot animal ~ plant                         #
##################################################
# 0. set up
	# 0.1. set working directory
		rm(list=ls()); gc()
		wd <- "~/Drive/COOPA_Analysis"
		data.files <- list.files(path = paste0(wd, "/Simulations/Data"))

# 1. plot difference ~ abiotic
	# 1.0. load and analyze data
		load(file = paste0(wd, "/Simulations/Data/", "diff.abiotic.Rdata"))
		sims <- dim(diff.abiotic)[1]
		n.resp <- dim(diff.abiotic)[2]
		n.pred <- dim(diff.abiotic)[3]
	# 1.1. plot and save .pdf
	# pdf(file = paste0(wd, "/Simulations/Figures/diff.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
	par(mfrow = c(n.resp, n.pred), mar = c(2, 2, 1, 0), oma = c(1, 0, 1, 1))
	for(i in 1:n.resp) {
		for (j in 1:n.pred) {
			obs.cor <- diff.abiotic[1, i, j]
			dens <- density(diff.abiotic[2:sims, i, j])
			x.min <- min(obs.cor, dens$x)
			x.max <- max(obs.cor, dens$x)
			plot(dens, xlim = c(x.min, x.max), ann = F, yaxt = "n")
			arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = max(dens$y)*0.9, col = "red", lwd = 2, length = 0)
			if (j == 1) {mtext(side = 2, text = dimnames(diff.abiotic)[[2]][i], cex = 0.75)}
			if (i == 1) {mtext(side = 3, text = dimnames(diff.abiotic)[[3]][j], cex = 0.75)}
			p.val <- round(
			length(which(abs(diff.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1)
			, digits = 2)
			if (p.val != 0) {
				legend("top", legend = bquote(hat("p")==.(p.val)), adj = c(0,.2), bty = "n")	
			} else {
				legend("top", legend = bquote(hat("p")==0), adj = c(0,0.2), bty = "n")
			}
		} # end of predictors, j
	} # end of responses, i
	# dev.off()

# 2. plot animal ~ abiotic
	# 2.0. load and analyze data
		load(file = paste0(wd, "/Simulations/Data/", "animal.abiotic.Rdata"))
		sims <- dim(animal.abiotic)[1]
		n.resp <- dim(animal.abiotic)[2]
		n.pred <- dim(animal.abiotic)[3]
	# 2.1. plot and save .pdf
		pdf(file = paste0(wd, "/Simulations/Figures/animal.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(2, 2, 1, 0), oma = c(1, 0, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- animal.abiotic[1, i, j]
				dens <- density(animal.abiotic[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				plot(dens, xlim = c(x.min, x.max), ann = F, yaxt = "n")
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = max(dens$y)*0.9, col = "red", lwd = 2, length = 0)
				if (j == 1) {mtext(side = 2, text = dimnames(animal.abiotic)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(animal.abiotic)[[3]][j], cex = 0.75)}
				p.val <- round(length(which(abs(animal.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				if (p.val != 0) {
					legend("top", legend = bquote(hat("p")==.(x.max)), adj = c(0,.2), bty = "n")	
				} else {
					legend("top", legend = bquote(hat("p")==0), adj = c(0,0.2), bty = "n")
				}		
			} # end of predictors, j
		} # end of responses, i
		dev.off()

# 3. plot plant ~ abiotic
	# 3.0. load and analyze data
		load(file = paste0(wd, "/Simulations/Data/", "plant.abiotic.Rdata"))
		sims <- dim(plant.abiotic)[1]
		n.resp <- dim(plant.abiotic)[2]
		n.pred <- dim(plant.abiotic)[3]
	# 3.1. plot and save .pdf
		pdf(file = paste0(wd, "/Simulations/Figures/plant.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(2, 2, 1, 0), oma = c(1, 0, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- plant.abiotic[1, i, j]
				dens <- density(plant.abiotic[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				plot(dens, xlim = c(x.min, x.max), ann = F, yaxt = "n")
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = max(dens$y)*0.9, col = "red", lwd = 2, length = 0)
				if (j == 1) {mtext(side = 2, text = dimnames(plant.abiotic)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(plant.abiotic)[[3]][j], cex = 0.75)}
				p.val <- round(length(which(abs(plant.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				if (p.val != 0) {
					legend("top", legend = bquote(hat("p")==.(x.max)), adj = c(0,.2), bty = "n")	
				} else {
					legend("top", legend = bquote(hat("p")==0), adj = c(0,0.2), bty = "n")
				}		
			} # end of predictors, j
		} # end of responses, i
		dev.off()

# 4. plot animal ~ plant
	# 4.0. load and analyze data
		load(file = paste0(wd, "/Simulations/Data/", "animal.plant.Rdata"))
		sims <- dim(animal.plant)[1]
		n.resp <- dim(animal.plant)[2]
		n.pred <- dim(animal.plant)[3]
	# 4.1. plot and save .pdf
		pdf(file = paste0(wd, "/Simulations/Figures/animal.plant.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(2, 2, 1, 0), oma = c(1, 0, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- animal.plant[1, i, j]
				dens <- density(animal.plant[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				plot(dens, xlim = c(x.min, x.max), ann = F, yaxt = "n")
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = max(dens$y)*0.9, col = "red", lwd = 2, length = 0)
				if (j == 1) {mtext(side = 2, text = dimnames(animal.plant)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(animal.plant)[[3]][j], cex = 0.75)}
				p.val <- round(length(which(abs(animal.plant[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				if (p.val != 0) {
					legend("top", legend = bquote(hat("p")==.(x.max)), adj = c(0,.2), bty = "n")	
				} else {
					legend("top", legend = bquote(hat("p")==0), adj = c(0,0.2), bty = "n")
				}		
			} # end of predictors, j
		} # end of responses, i
		dev.off()