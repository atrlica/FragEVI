library(raster)
library(data.table)


bos.1m <- stack("")



## for loss projections, could take the median collection in each pixel and 
## throw dice for each member given its size or age-based mortality expectation
## do this 20 times to get the likely mortality for each set
## each year there is a mortality replace with 5cm tree and let grow at expected pace
## and subject it to mortality risk every year as you go.
## You could take a collection of the succesful simulation collections (say within 5% of the target biomass)
## select 100 at random and subject each to this 20 year mortality test and find out what the final
## biomass + productivity will be.

### for gain projection, you could (while allowing the mortalities above to take shape)
### ID places that are 1) pervious + 2) not canopied
### and then somehow identify contiguous patches of ground that will allow the average
### local tree to get to say 30cm -- could use street tree alometrics to figure out what
## a reasonable canopy area this would represent.
pi*(10^2)