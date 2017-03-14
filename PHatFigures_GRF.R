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
		data.wd <- paste0(wd, "/Simulations/Data/GRFData")
		figures.wd <- paste0(wd, "/Simulations/Figures/GRFFigs")
		data.files <- list.files(path = data.wd, pattern = ".Rdata")

# 1. plot difference ~ abiotic
	# 1.0. load and analyze data
		load(file = paste0(data.wd, "/diff.abiotic.Rdata"))
		sims <- dim(diff.abiotic)[1]
		n.resp <- dim(diff.abiotic)[2]
		n.pred <- dim(diff.abiotic)[3]
	# 1.1. plot and save .pdf
	pdf(file = paste0(figures.wd, "/diff.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(1, 1, 0, 0), oma = c(3, 1, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- diff.abiotic[1, i, j]
				dens <- density(diff.abiotic[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				y.min <- min(obs.cor, dens$y)
				y.max <- max(obs.cor, dens$y)
				x.vals <- dens$x
				y.vals <- dens$y
				x.extreme <- x.vals[x.vals > obs.cor]
				y.extreme <- y.vals[x.vals > obs.cor]
				x.extreme.low <- x.vals[x.vals < obs.cor]
				y.extreme.low <- y.vals[x.vals < obs.cor]
				plot(NA, xlim = c(-1, 1), ylim = c(0, (y.max*1.25)), ann = F, yaxt = "n", xaxt = "n")
				axis(1, labels = F)
				if (j == 1) {mtext(side = 2, text = dimnames(diff.abiotic)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(diff.abiotic)[[3]][j], cex = 0.75)}
				if (i == n.resp) {axis(1, tick = F, labels = T)}
				p.val <- round(length(which(abs(diff.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				larger.cors <- c(length(which(diff.abiotic[2:sims, i, j] < obs.cor)), length(which(diff.abiotic[2:sims, i, j] > obs.cor)))
				if (larger.cors[2] > larger.cors[1]) {
					p.val <- larger.cors[1]/sims
				} else {
					p.val <- larger.cors[2]/sims
				}
				colorramp <- colorRampPalette(c("blue", "red"))
				if (larger.cors[1] > larger.cors[2]) {
					poly.col <- colorramp(sims)[larger.cors[1]]
					polygon(x = c(x.extreme,rev(x.extreme)), y = c(rep(0,length(x.extreme)), rev(y.extreme)), col = poly.col, border = F)
				} else {
					poly.col <- colorramp(sims)[larger.cors[2]]
					polygon(x = c(x.extreme.low,rev(x.extreme.low)), y = c(rep(0,length(x.extreme.low)), rev(y.extreme.low)), col = poly.col, border = F)
				}
				lines(x = x.vals, y = y.vals)
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = y.max, col = "black", lwd = 1.5, length = 0)				
				if (p.val <= 0.05) {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val,2))), adj = 0, col = "red")
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0, col = "red")		
					} else {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val, 2))), adj = 0)
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0)		
				}
			} # end of predictors, j
		} # end of responses, i
		mtext(outer = T, text = bquote("Spearman's correlation coefficient,"~rho), line = 2, side = 1)
	dev.off()

# 2. plot animal ~ abiotic
	# 2.0. load and analyze data
		load(file = paste0(data.wd, "/animal.abiotic.Rdata"))
		sims <- dim(animal.abiotic)[1]
		n.resp <- dim(animal.abiotic)[2]
		n.pred <- dim(animal.abiotic)[3]
	# 2.1. plot and save .pdf
		pdf(file = paste0(figures.wd, "/animal.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(1, 1, 0, 0), oma = c(3, 1, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- animal.abiotic[1, i, j]
				dens <- density(animal.abiotic[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				y.min <- min(obs.cor, dens$y)
				y.max <- max(obs.cor, dens$y)
				x.vals <- dens$x
				y.vals <- dens$y
				x.extreme <- x.vals[x.vals > obs.cor]
				y.extreme <- y.vals[x.vals > obs.cor]
				x.extreme.low <- x.vals[x.vals < obs.cor]
				y.extreme.low <- y.vals[x.vals < obs.cor]
				plot(NA, xlim = c(-1, 1), ylim = c(0, (y.max*1.25)), ann = F, yaxt = "n", xaxt = "n")
				axis(1, labels = F)
				if (j == 1) {mtext(side = 2, text = dimnames(animal.abiotic)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(animal.abiotic)[[3]][j], cex = 0.75)}
				if (i == n.resp) {axis(1, tick = F, labels = T)}
				p.val <- round(length(which(abs(animal.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				larger.cors <- c(length(which(animal.abiotic[2:sims, i, j] < obs.cor)), length(which(animal.abiotic[2:sims, i, j] > obs.cor)))
				if (larger.cors[2] > larger.cors[1]) {
					p.val <- larger.cors[1]/sims
				} else {
					p.val <- larger.cors[2]/sims
				}
				colorramp <- colorRampPalette(c("blue", "red"))
				if (larger.cors[1] > larger.cors[2]) {
					poly.col <- colorramp(sims)[larger.cors[1]]
					polygon(x = c(x.extreme,rev(x.extreme)), y = c(rep(0,length(x.extreme)), rev(y.extreme)), col = poly.col, border = F)
				} else {
					poly.col <- colorramp(sims)[larger.cors[2]]
					polygon(x = c(x.extreme.low,rev(x.extreme.low)), y = c(rep(0,length(x.extreme.low)), rev(y.extreme.low)), col = poly.col, border = F)
				}
				lines(x = x.vals, y = y.vals)
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = y.max, col = "black", lwd = 1.5, length = 0)				
				if (p.val <= 0.05) {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val,2))), adj = 0, col = "red")
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0, col = "red")		
					} else {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val, 2))), adj = 0)
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0)		
				}
				} # end of predictors, j
		} # end of responses, i
		mtext(outer = T, text = bquote("Spearman's correlation coefficient,"~rho), line = 2, side = 1)
		dev.off()

# 3. plot plant ~ abiotic
	# 3.0. load and analyze data
		load(file = paste0(data.wd, "/plant.abiotic.Rdata"))
		sims <- dim(plant.abiotic)[1]
		n.resp <- dim(plant.abiotic)[2]
		n.pred <- dim(plant.abiotic)[3]
	# 3.1. plot and save .pdf
		pdf(file = paste0(figures.wd, "/plant.abiotic.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(1, 1, 0, 0), oma = c(3, 1, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- plant.abiotic[1, i, j]
				dens <- density(plant.abiotic[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				y.min <- min(obs.cor, dens$y)
				y.max <- max(obs.cor, dens$y)
				x.vals <- dens$x
				y.vals <- dens$y
				x.extreme <- x.vals[x.vals > obs.cor]
				y.extreme <- y.vals[x.vals > obs.cor]
				x.extreme.low <- x.vals[x.vals < obs.cor]
				y.extreme.low <- y.vals[x.vals < obs.cor]
				plot(NA, xlim = c(-1, 1), ylim = c(0, (y.max*1.25)), ann = F, yaxt = "n", xaxt = "n")
				axis(1, labels = F)
				if (j == 1) {mtext(side = 2, text = dimnames(plant.abiotic)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(plant.abiotic)[[3]][j], cex = 0.75)}
				if (i == n.resp) {axis(1, tick = F, labels = T)}
				p.val <- round(length(which(abs(plant.abiotic[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				larger.cors <- c(length(which(plant.abiotic[2:sims, i, j] < obs.cor)), length(which(plant.abiotic[2:sims, i, j] > obs.cor)))
				if (larger.cors[2] > larger.cors[1]) {
					p.val <- larger.cors[1]/sims
				} else {
					p.val <- larger.cors[2]/sims
				}
				colorramp <- colorRampPalette(c("blue", "red"))
				if (larger.cors[1] > larger.cors[2]) {
					poly.col <- colorramp(sims)[larger.cors[1]]
					polygon(x = c(x.extreme,rev(x.extreme)), y = c(rep(0,length(x.extreme)), rev(y.extreme)), col = poly.col, border = F)
				} else {
					poly.col <- colorramp(sims)[larger.cors[2]]
					polygon(x = c(x.extreme.low,rev(x.extreme.low)), y = c(rep(0,length(x.extreme.low)), rev(y.extreme.low)), col = poly.col, border = F)
				}
				lines(x = x.vals, y = y.vals)
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = y.max, col = "black", lwd = 1.5, length = 0)				
				if (p.val <= 0.05) {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val,2))), adj = 0, col = "red")
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0, col = "red")		
					} else {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val, 2))), adj = 0)
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0)		
				}	
			} # end of predictors, j
		} # end of responses, i
		mtext(outer = T, text = bquote("Spearman's correlation coefficient,"~rho), line = 2, side = 1)
		dev.off()

# 4. plot animal ~ plant
	# 4.0. load and analyze data
		load(file = paste0(data.wd, "/animal.plant.Rdata"))
		sims <- dim(animal.plant)[1]
		n.resp <- dim(animal.plant)[2]
		n.pred <- dim(animal.plant)[3]
	# 4.1. plot and save .pdf
		pdf(file = paste0(figures.wd, "/animal.plant.pdf"), width = n.pred, height = n.resp, pointsize = 8)
		par(mfrow = c(n.resp, n.pred), mar = c(1, 1, 0, 0), oma = c(3, 1, 1, 1))
		for(i in 1:n.resp) {
			for (j in 1:n.pred) {
				obs.cor <- animal.plant[1, i, j]
				dens <- density(animal.plant[2:sims, i, j])
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				x.min <- min(obs.cor, dens$x)
				x.max <- max(obs.cor, dens$x)
				y.min <- min(obs.cor, dens$y)
				y.max <- max(obs.cor, dens$y)
				x.vals <- dens$x
				y.vals <- dens$y
				x.extreme <- x.vals[x.vals > obs.cor]
				y.extreme <- y.vals[x.vals > obs.cor]
				x.extreme.low <- x.vals[x.vals < obs.cor]
				y.extreme.low <- y.vals[x.vals < obs.cor]
				plot(NA, xlim = c(-1, 1), ylim = c(0, (y.max*1.25)), ann = F, yaxt = "n", xaxt = "n")
				axis(1, labels = F)
				if (j == 1) {mtext(side = 2, text = dimnames(animal.plant)[[2]][i], cex = 0.75)}
				if (i == 1) {mtext(side = 3, text = dimnames(animal.plant)[[3]][j], cex = 0.75)}
				if (i == n.resp) {axis(1, tick = F, labels = T)}
				p.val <- round(length(which(abs(animal.plant[2:sims, i, j]) > abs(obs.cor)) == T)/(sims-1), digits = 2)
				larger.cors <- c(length(which(animal.plant[2:sims, i, j] < obs.cor)), length(which(animal.plant[2:sims, i, j] > obs.cor)))
				if (larger.cors[2] > larger.cors[1]) {
					p.val <- larger.cors[1]/sims
				} else {
					p.val <- larger.cors[2]/sims
				}
				colorramp <- colorRampPalette(c("blue", "red"))
				if (larger.cors[1] > larger.cors[2]) {
					poly.col <- colorramp(sims)[larger.cors[1]]
					polygon(x = c(x.extreme,rev(x.extreme)), y = c(rep(0,length(x.extreme)), rev(y.extreme)), col = poly.col, border = F)
				} else {
					poly.col <- colorramp(sims)[larger.cors[2]]
					polygon(x = c(x.extreme.low,rev(x.extreme.low)), y = c(rep(0,length(x.extreme.low)), rev(y.extreme.low)), col = poly.col, border = F)
				}
				lines(x = x.vals, y = y.vals)
				arrows(x0 = obs.cor, x1 = obs.cor, y0 = 0, y1 = y.max, col = "black", lwd = 1.5, length = 0)				
				if (p.val <= 0.05) {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val,2))), adj = 0, col = "red")
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0, col = "red")		
					} else {
					text(x = -1, y = y.max*1.15, labels = bquote(hat("p")==.(round(p.val, 2))), adj = 0)
					text(x = -1, y = y.max*0.95, labels = bquote(rho^"*"==.(round(obs.cor, 2))), adj = 0)		
				}	
			} # end of predictors, j
		} # end of responses, i
		mtext(outer = T, text = bquote("Spearman's correlation coefficient,"~rho), line = 2, side = 1)
		dev.off()