# ToxicoDI

An R package for building ToxicoDB database tables from a list of ToxicoSet objects.

## Usage

```{r}
library(ToxicoDI)

# Get individual tables from each TSet and write them to procdata
extractAllTables(path='tSets',
    pattern='EM.*rds|.*ldh.*rds|.*drugMatrix.*rds', 
    outDir='procdata')

# Perform joins on  individual PSet tables to get the final database tables.
buildAllTables()

# Results will be written to the 'latest' directory.

```