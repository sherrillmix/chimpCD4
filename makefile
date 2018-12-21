all: out/glm.pdf out/cpzFits.pdf out/monkeyFits.pdf out/mbFits.pdf

out/glm.pdf:
	R CMD BATCH analyzeSites.R

out/cpzFits.pdf:
	R CMD BATCH analyzeCpz.R

out/monkeyFits.pdf:
	R CMD BATCH analyzeMonkey.R

out/mbFits.pdf:
	R CMD BATCH analyzeMB897.R

