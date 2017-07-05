# GRB_interactions

The analysis runs in two steps.

### 1. Download and format data 

Note, that you have to copy the GRB, CNE, and
target/bystander gene data (GRB_Interactions.zip) into the folder
`data/GRB_Interactions/` because they are not part of this repository.

Than execute the following:
``` 
cd data 
sh download.sh 
```


### 2. Run analysis 

``` 
Rscript R/GRB_interactions.R 
```