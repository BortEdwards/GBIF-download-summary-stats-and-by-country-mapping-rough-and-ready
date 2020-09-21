##### Quick and dirty code to download data from GBIF, prune, produce basic by-country summary statistics, and map results.
##### Requires list of countries with areas (from Wikipedia: https://en.wikipedia.org/wiki/List_of_countries_and_dependencies_by_area)
##### Mapping uses base R functions or plotly and could/should be upgraded to ggplot2
##### Requires user to have a GBIF account
##### https://docs.ropensci.org/rgbif/articles/downloads.html - tutorial for downloading from GBIF
##### WARNING!! There are no data cleaning or taxonomy harmonization steps between GBIF and finished product!!!


library(rgbif)
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_
library(tidyr)
library(countrycode) # to convert between names and ISC country codes

user <- "***"
pwd <- "***"
email <- "***"
# ideally login details should be saved in R environment (safer)

# key <- name_backbone(name='Puma concolor')$speciesKey  # find key for desired taxon
# key <- 3065  # can call on list of multiple taxa if required

# sends request to GBIF to start compiling query for download
GBIF_DL <- occ_download(
  pred_in("taxonKey", 3065), # if key not designated above, nominate here
  pred("hasCoordinate", TRUE), # define requirements for data here
  user = user, pwd = pwd, email = email
  )

# check on status of query
occ_download_meta(GBIF_DL)

#occ_download_cancel(key="***", user = "***", pwd = "***") # if need to cancel GBIF query

# When query successfully complete, it can be retrieved either by entering link in browser, or directly via R (changing ID number):

dat <- occ_download_get(key="***", overwrite=TRUE) # Key # retrieved above
dat2 <- occ_download_import(dat)

as.data.frame(colnames(dat2)) #lists column names

subset.dat <- data.frame(dat2$countryCode, dat2$species, dat2$speciesKey) #extract columns of interest (this is only done to reduce dataframe size for limited computation space)
# dat2$organismID, dat2$organismName, dat2$identificationID, dat2$taxonID*, dat2$acceptedNameUsage, dat2$previousIdentifications, dat2$nameAccordingTo # removed as mostly empty

-------
  
subset.dat2 <- subset.dat %>% # rename columns
  rename(country = dat2.countryCode,
         species = dat2.species)

subset.dat2 <- filter(subset.dat2, dat2.speciesKey != "NA") # remove rows where species key is NA (these are generic only IDs)

subset.dat3 <- subset(subset.dat2, select = -dat2.speciesKey) # remove species key column (this is only done to reduce dataframe size for limited computation space)

write.csv(subset.dat3,'SpeciesByCountry.csv')

numcntry.records <- count(subset.dat3, country) # number of records per country
numcntry.records <- numcntry.records %>% # rename new column
  rename(
    RecordsPerCountry = n
  )

numcntry.species <- subset.dat3 [!duplicated(subset.dat3[,c("country", "species")]),] # reduce to one record per country
numcntry.species.tally <- numcntry.species %>% group_by(country) %>% tally() # number of species for each country

numcntry.species.tally <- numcntry.species.tally %>% # rename columns
              rename(
              TotalSpecies = n
              )

numcntry.unique <- numcntry.species %>% # filter out species that occur in more than one country
  group_by(species) %>% 
  filter(n()==1)

#write.csv(subset.dat4,'UniqueSpeciesByCountry.csv')

numcntry.endemic <- numcntry.unique %>% group_by(country) %>% tally() # tally number of unique species for each country

numcntry.endemic <- numcntry.endemic %>% # rename columns
              rename(
              EndemicSpecies = n
              )

Spp.Totals <- numcntry.species.tally %>% full_join(numcntry.endemic) %>% # merge total records, total species, and endemic species dataframes - introducing NAs
full_join(numcntry.records) %>%
  mutate_if(is.numeric,coalesce,0)  # replace NA with zeros for numeric columns

Spp.Totals <- within(Spp.Totals, ProportionEndemicofTotalSpp <- EndemicSpecies/TotalSpecies) # add column with proportion of species that are endemic

Spp.Totals$iso3c   <- countrycode(Spp.Totals$country, origin = "iso2c", destination = "iso3c") # add ISO3 codes

write.csv(Spp.Totals,'TotalandEndemicSpp.csv')

---------
  
##### Add country size and calculate species proportions

df1 <- read.csv("CountrySizes2.csv") # import csv with countries and associated areas
df1 <- df1 %>% drop_na


df1$iso2c   <- countrycode(df1$country, origin = "country.name", destination = "iso2c") # Convert country names to ISO2 and ISO3 codes
df1 <- select(df1, km, iso2c) # subset dataframe
df1 <- rename(df1, c("iso2c"="country")) # standardize column names
df1 <- df1 %>% drop_na

TotalandEndemicSppCalcs <- Spp.Totals %>% full_join(df1) # merge country dataframe values with main dataframe

TotalandEndemicSppCalcs <- within(TotalandEndemicSppCalcs, SppPerKm <- TotalSpecies/km) # add column with species per km squares area of country

TotalandEndemicSppCalcs <- within(TotalandEndemicSppCalcs, EndemicSppPerKm <- EndemicSpecies/km) # add column with number of endemic species per kms for each country

TotalandEndemicSppCalcs <- within(TotalandEndemicSppCalcs, RecordsPerKm <- RecordsPerCountry/km) # add column with number of records per km2 for each country

TotalandEndemicSppCalcs <- within(TotalandEndemicSppCalcs, RecordsPerSpecies <- RecordsPerCountry/TotalSpecies) # add column with avg number of records per species for each country

write.csv(TotalandEndemicSppCalcs,'TotalandEndemicSppWithCalcs.csv')

---------
  
###### PLOT VALUES PER COUNTRY ON WORLD MAP
###### https://stackoverflow.com/questions/24136868/plot-map-with-values-for-countries-as-color-in-r

library(ggmap)
library(RColorBrewer)
library(maptools)
library(plyr)

#THIS IS WHERE YOU NOMINATE THE COLUMN TO PLOT ON THE MAP BY CHANGING ITS NAME TO "VALUE"
ddf <- rename(TotalandEndemicSppCalcs, c("EndemicSpecies"="value")) # rename focal variable column

# the below removes rows with missing data that doesn't play nice with mapping
ddf <- ddf[!(ddf$country == ""), ] # remove rows with no country code
ddf <- ddf[!(ddf$country == "ZZ"),] # remove rows with GBIF ambiguity code "ZZ"
ddf <- ddf[!(ddf$country == "BQ"),] # remove rows with "BQ" (Caribbean "special municipalities")
ddf <- ddf[!(ddf$country == "NA"),] # remove rows with "NA" missing data in "country" field
ddf <- ddf[!(ddf$country == "AQ"),] # remove Antarctica
ddf <- ddf[!(ddf$country == "CW"),] # remove Curaçao as missing from map data
ddf <- ddf[!(ddf$country == "SX"),] # remove Sint Maarten as missing from map data
#ddf <- ddf[!(ddf$km == "NA"),] #remove rows with "NA" missing data in "km" field only necessary for calculations/maps involving country areas
ddf <- na.omit(ddf) #remove empty rows
ddf <- as.data.frame(ddf)

data(wrld_simpl) # load base world map

plotSppPerCountry <- function() {
  
  wrld <- wrld_simpl[wrld_simpl@data$UN!="10",] # remove Antarctica

  pal <- colorRampPalette(brewer.pal(9, 'Reds'))(75) # noninate color scheme and number of breaks along spectrum
  pal <- pal[with(ddf, findInterval(value, sort(unique(value))))]
  col <- rep(grey(0.8), length(wrld_simpl@data$ISO2))
  col[match(ddf$country, wrld@data$ISO2)] <- pal
  plot(wrld, col = col)
 
}

plotSppPerCountry()


##### Alternative mapping package with built-in scale etc
##### https://plotly.com/r/choropleth-maps/

library(plotly)

l <- list(color = toRGB("grey"), width = 0.5) # light grey boundaries

g <- list(
  showframe = FALSE,
  showcoastlines = FALSE,
  projection = list(type = 'Mercator') # specify map projection/options
)

fig <- plot_geo(ddf)
fig <- fig %>% add_trace(
  z = ~value, color = ~value, colors = 'Reds',
  text = ~iso3c, locations = ~iso3c, marker = list(line = l)
)
fig <- fig %>% colorbar(title = 'LegendTitle', tickprefix = 'Spp')
fig <- fig %>% layout(
  title = 'TopTitle',
  geo = g
)

fig

