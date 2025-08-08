# Assignment: Muhammad Ekramah - Class Ib
#

# 1. Project Setup and Folder Organization

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# 2. Data Loading and Initial Inspection

csv_file_path <- "/home/ekramah/AI and Omics Research Internship (2025)/Module I: Getting Started with R (Class Ib)/Module_I/patient_info.csv"

# Load the data using the absolute path
patient_data <- read.csv(csv_file_path)


str(patient_data)

View(patient_data)

# 3. Data Cleaning and Type Conversion
# Convert character columns to factors for proper categorical analysis.

# Convert 'diagnosis' to a factor with specific levels
patient_data$diagnosis_fact <- as.factor(patient_data$diagnosis)
patient_data$diagnosis_fact <- factor(patient_data$diagnosis_fact,
                                      levels = c("Normal", "Cancer"))
# Check the levels
levels(patient_data$diagnosis_fact)

# Convert 'gender' to a factor
patient_data$gender_fact <- as.factor(patient_data$gender)
patient_data$gender_fact <- factor(patient_data$gender_fact,
                                   levels = c("Male", "Female"))
# Check the levels
levels(patient_data$gender_fact)

# Create a new binary numeric variable for gender (Male = 1, Female = 0)
patient_data$gender_num <- ifelse(patient_data$gender_fact == "Male", 1, 0)
class(patient_data$gender_num)

# Convert 'smoker' to a factor
patient_data$smoker_fact <- as.factor(patient_data$smoker)
levels(patient_data$smoker_fact)

# Create a new binary variable for smoking status (Yes = 1, No = 0)
patient_data$smoking_status <- ifelse(patient_data$smoker_fact == "Yes", 1, 0)
class(patient_data$smoking_status)

# Final check of the data structure after conversions
str(patient_data)

# 4. Saving the Outputs
# Save the cleaned data frame to a new CSV file in the 'clean_data' folder.

write.csv(patient_data,
          file = file.path("clean_data", "patient_info_clean.csv"),
          row.names = FALSE)

