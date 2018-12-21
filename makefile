all:
	R CMD BATCH analyzeSites.R
	R CMD BATCH analyzeCpz.R
	R CMD BATCH analyzeMB897.R
	R CMD BATCH analyzeMonkey.R
