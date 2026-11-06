# Fill design_template.csv with your four justified alternative settings.
# Each of exactly two factors has two levels; keep the other two factors fixed.
# Run from the student packet folder, after renaming your filled table design.csv.
source("hw2_tools.R")
plan <- read.csv("design.csv",stringsAsFactors=FALSE)
# Set B and seed explicitly; retain this call in your submitted analysis script.
sim <- run_plan(plan,B=1000,seed=81082026,output_dir="results/simulation")
# Add your plots and Task 2 analysis. This script is a starting point, not a submission.
