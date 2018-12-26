all: out/glm.pdf out/cpzFits.pdf out/monkeyFits.pdf out/mbFits.pdf

out/glm.pdf:
	R CMD BATCH --no-save --no-restore analyzeSites.R

out/cpzFits.pdf:
	R CMD BATCH --no-save --no-restore analyzeCpz.R

out/monkeyFits.pdf:
	R CMD BATCH --no-save --no-restore analyzeMonkey.R

out/mbFits.pdf:
	R CMD BATCH --no-save --no-restore analyzeMB897.R

