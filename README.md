# coinfection_heterogeneity
Corresponds with the article titled "Heterogeneity in viral infection increases the rate of deleterious mutation accumulation"


To reproduce the data from each figure, simply run the appropriate "params_" file.
To reproduce the figures from the processed date, run the appropriate "Visualize " R Markdown file.

Note that to recapitulate the models where mutations are treated as either recessive or dominant, 
one needs to replace the files in this main directory with their counterparts in the respective 
recessive or dominant directories.

The "main_coinfection_script" executes the "run_N_script" with the appropriate parameterization and 
model assumptions. It also also controls the number of replicates that are performed. The run_N_script
iterates through each viral generation and saves the data for that replicate. The "summarize_data_NeV"
script reformats the data by taking the mean number of mutations accumulated per virion per generation 
per replicate. The summary output files are then used in the Visualization R files.

Please contact me if you have questions or want advice on how to modify the code to change the 
biological context or assumptions.
