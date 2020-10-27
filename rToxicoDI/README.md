# ToxicoDB
ToxicoDB web app

## Database Tables

This branch holds the scripts used to regenerate the ToxicoDB database tables from ToxicoSet R objects using the ToxicoGx R package.

### Instructions

1. Download the tSets you wish to use into the tSets directory
2. Open R or Rstudio
3. Change into the scripts directory: `setwd(scripts)`
4. Run `O_runTableGenerationPipeline.R`
5. Final csv files of tables are tagged `_latest.csv`
6. Upload to MySQL using your migration scripts
