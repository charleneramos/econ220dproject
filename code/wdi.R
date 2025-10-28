# =============================================================================#
#   Setup R for World Bank World Development Indicators (WDI)                  # 
#   https://databank.worldbank.org/source/world-development-indicators         #
# =============================================================================#

cat("\014")
rm(list = ls())

install.packages("WDI") # https://github.com/vincentarelbundock/WDI
install.packages("wbstats") # https://econandrew.github.io/wdi-api-R-gettingstarted/using-r-to-access-the-world-bank-api.html
install.packages("magick")
install.packages("vars")

library(WDI)
library(wbstats) 
library(magick)
library(vars)

# =============================================================================#
#   II(a) Search for relevant macroeconomic variables                                # 
# =============================================================================#

# Get a list of countries and identify argument entry for Philippines
countries <- WDI_data$country # PH

# GDP (current US$)
gdp <- WDIsearch(string = 'gdp', field = 'name', cache = NULL)
ph_gdp <- WDI(country = 'PH', indicator = 'NY.GDP.MKTP.CD')
sum(!is.na(ph_gdp$NY.GDP.MKTP.CD)) # 65

# Unemployment, total (% of total labor force) (modeled ILO estimate)
unrate <- WDIsearch(string = 'unemployment', field = 'name', cache = NULL)
ph_unrate <- WDI(country = 'PH', indicator = 'SL.UEM.TOTL.ZS')
sum(!is.na(ph_unrate$SL.UEM.TOTL.ZS)) # 34

# Inflation, consumer prices (annual %)
cpi <- WDIsearch(string = 'inflation', field = 'name', cache = NULL)
ph_cpi <- WDI(country = 'PH', indicator = 'FP.CPI.TOTL.ZG')
sum(!is.na(ph_cpi$FP.CPI.TOTL.ZG)) # 14

# Population, total
pop <- WDIsearch(string = 'population', field = 'name', cache = NULL)
ph_pop <- WDI(country = 'PH', indicator = 'SP.POP.TOTL')
sum(!is.na(ph_pop$SP.POP.TOTL)) # 65

# Labor force, total
lf <- WDIsearch(string = 'labor force', field = 'name', cache = NULL)
ph_lf <- WDI(country = 'PH', indicator = 'SL.TLF.TOTL.IN')
sum(!is.na(ph_lf$SL.TLF.TOTL.IN)) # 35

# Foreign direct investment, net inflows (BoP, current US$)
fdi <- WDIsearch(string = 'foreign direct investment', field = 'name', cache = NULL)
ph_fdi <- WDI(country = 'PH', indicator = 'BX.KLT.DINV.CD.WD')
sum(!is.na(ph_fdi$BX.KLT.DINV.CD.WD)) # 55 

# Net official development assistance received (current US$)
oda <- WDIsearch(string = 'Net official development assistance received', field = 'name', cache = NULL)
ph_oda <- WDI(country = 'PH', indicator = 'DT.ODA.ODAT.CD')
sum(!is.na(ph_oda$DT.ODA.ODAT.CD)) # 64 

# Personal remittances, received (current US$)
remit <- WDIsearch(string = 'personal remittances', field = 'name', cache = NULL)
ph_remit <- WDI(country = 'PH', indicator = 'BX.TRF.PWKR.CD.DT')
sum(!is.na(ph_remit$BX.TRF.PWKR.CD.DT)) # 48

rm(list = c('cpi', 'fdi', 'gdp', 'lf', 'oda', 'pop', 'remit', 'unrate'))
rm(list = c('ph_unrate','ph_cpi','ph_lf'))

# =============================================================================#
#   II(b) Create merged dataset                                                      # 
# =============================================================================#

# Remove NA values
ph_fdi <- na.omit(ph_fdi)
ph_gdp <- na.omit(ph_gdp)
ph_oda <- na.omit(ph_oda)
ph_pop <- na.omit(ph_pop)
ph_remit <- na.omit(ph_remit)

# Reorder by ascending year
ph_fdi <- ph_fdi[order(ph_fdi$year), ]
ph_gdp <- ph_gdp[order(ph_gdp$year), ]
ph_oda <- ph_oda[order(ph_oda$year), ]
ph_pop <- ph_pop[order(ph_pop$year), ]
ph_remit <- ph_remit[order(ph_remit$year), ]

# Drop irrelevant columns
ph_fdi <- ph_fdi[, !(names(ph_fdi) %in% c('country','iso3c'))]
ph_gdp <- ph_gdp[, !(names(ph_gdp) %in% c('country','iso3c'))]
ph_oda <- ph_oda[, !(names(ph_oda) %in% c('country','iso3c'))]
ph_pop <- ph_pop[, !(names(ph_pop) %in% c('country','iso3c'))]
ph_remit <- ph_remit[, !(names(ph_remit) %in% c('country','iso3c'))]

# Rename columns
names(ph_fdi) <- c('country','year','fdi')
names(ph_gdp) <- c('country','year','gdp')
names(ph_oda) <- c('country','year','oda')
names(ph_pop) <- c('country','year','pop')
names(ph_remit) <- c('country','year','remit')

# Merge columns
db <- Reduce(function(x, y) merge(x, y, by = c('country', 'year')), 
             list(ph_fdi, ph_gdp, ph_oda, ph_pop, ph_remit))


db <- db[, c('year','pop','gdp','fdi','oda','remit')]

# Short to long dataset
db_long <- reshape(db,
                    varying = list(c('pop','gdp','fdi','oda','remit')),
                    v.names = 'value',
                    timevar = 'variable',
                    times = c('pop','gdp','fdi','oda','remit'),
                    direction = 'long')

db_long <- db_long[, !(names(db_long) == "id")]
rownames(db_long) <- NULL

rm('ph_fdi','ph_gdp','ph_oda','ph_pop','ph_remit')

# =============================================================================#
#   II(c) Plot assembled data                                                        # 
# =============================================================================#

# Principal component analysis
pca <- prcomp(db[, c('pop','gdp','fdi','oda','remit')],
              center = TRUE, scale. = TRUE)
summary(pca)

# Population, total
plot_pop <- db_long[db_long$variable == 'pop', ]
plot_pop$value_mil <- plot_pop$value / 1e6 
plot_pop$value_bil <- plot_pop$value / 1e9 

png('output/plot_pop.png', width = 800, height = 600, res = 120)

plot(plot_pop$year, plot_pop$value_mil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     main = 'Philippines Total Population, 1977-2023',
     xlab = 'Year', ylab = 'Population Levels (Millions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_pop <- image_read('output/plot_pop.png')
png_pop <- image_trim(png_pop)
image_write(png_pop, path = 'output/plot_pop.png')

# GDP (current US$)

plot_gdp <- db_long[db_long$variable == 'gdp', ]
plot_gdp$value_mil <- plot_gdp$value / 1e6 
plot_gdp$value_bil <- plot_gdp$value / 1e9 

png('output/plot_gdp.png', width = 800, height = 600, res = 120)

plot(plot_gdp$year, plot_gdp$value_bil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     main = 'Philippines GDP, 1977-2023',
     xlab = 'Year', ylab = 'Current US$ (Billions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_gdp <- image_read('output/plot_gdp.png')
png_gdp <- image_trim(png_gdp)
image_write(png_gdp, path = 'output/plot_gdp.png')

# Foreign direct investment, net inflows (BoP, current US$)

plot_fdi <- db_long[db_long$variable == 'fdi', ]
plot_fdi$value_mil <- plot_fdi$value / 1e6 
plot_fdi$value_bil <- plot_fdi$value / 1e9 

png('output/plot_fdi.png', width = 800, height = 600, res = 120)

plot(plot_fdi$year, plot_fdi$value_bil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     main = 'Philippines FDI, 1977-2023',
     xlab = 'Year', ylab = 'Current US$ (Billions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_fdi <- image_read('output/plot_fdi.png')
png_fdi <- image_trim(png_fdi)
image_write(png_fdi, path = 'output/plot_fdi.png')

# =============================================================================#
#   III VAR and IRF analysis                                                   # 
# =============================================================================#

# Prepare data for VAR



