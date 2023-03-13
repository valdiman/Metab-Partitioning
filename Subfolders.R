# Create subfolder within the current working directory
dir.create("Data")
dir.create("R")
dir.create("Output")
dir.create("Output/Data")
dir.create("Output/Data/csv")
dir.create("Output/Data/excel")
dir.create("Output/Plots")

# Move files
file.copy(from = "data.xlsx",
          to = "Data/data.xlsx")
# Delete files
file.remove(from = "data.xlsx")
