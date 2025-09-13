# 📘 R + AMR Package Cheat Sheet (Beginner to Advanced)

## 🧠 1. R BASICS FOR BIOINFORMATICS & AMR

### Variables and Data Types
```r
x <- 5                # Numeric
name <- "E. coli"     # Character
test <- TRUE          # Logical
```

### Vectors and Data Frames
```r
a <- c(1, 2, 3)              # Numeric vector
species <- c("E. coli", "S. aureus")
df <- data.frame(id = 1:2, species = species)  # Data frame
```

### Common Functions
```r
length(a)          # Number of elements
str(df)            # Structure of object
summary(df)        # Quick summary
class(x)           # Data type
```

### Conditional & Loops
```r
if (x > 0) print("Positive")
for (i in 1:5) print(i)
```

---

## 🧰 2. CORE `AMR` PACKAGE FUNCTIONS

### Load Package
```r
library(AMR)
```

### Microorganism Cleaning
```r
mo <- as.mo("E. coli")       # Get MO code
mo_fullname(mo)              # Full name: Escherichia coli
mo_gramstain(mo)             # Gram stain: Gram-negative
```

### Antibiotic Data Standardization
```r
as.ab("AMP")       # Recognizes Ampicillin
ab_name("AMP")     # Ampicillin
ab_group("AMP")    # Group: Beta-lactams
```

### S/I/R Interpretation
```r
as.rsi("S")                    # Converts to class 'rsi'
interpretive_reading(df)      # EUCAST rules to infer resistance
```

---

## 📊 3. RESISTANCE ANALYSIS

### Basic Resistance Rate
```r
resistance(df$AMP_rsi)        # % resistant
susceptibility(df$CIP_rsi)    # % susceptible
```

### Grouped Summary
```r
df %>% group_by(hospital) %>% 
  summarise(res = resistance(CIP_rsi))
```

---

## 🦠 4. MULTI-DRUG RESISTANCE (MDR/XDR/PDR)

```r
mdro(df, guideline = "EUCAST")  # Detect MDR based on EUCAST rules
```

### Example Output Categories:
- "Negative"
- "Positive"
- "Positive, unconfirmed"

---

## 📈 5. TREND & TIME-SERIES ANALYSIS

### Resistance Rate Over Time
```r
ts_df <- df %>% group_by(year) %>% 
  summarise(res = resistance(AMC_rsi))
```

### ARIMA Forecasting (Optional)
```r
library(forecast)
ts_obj <- ts(ts_df$res, start=2015, frequency=1)
fit <- auto.arima(ts_obj)
forecast(fit, h=3)
```

---

## 💾 6. EXPORT RESULTS
```r
write.csv(df, "results.csv")         # Export data
ggsave("plot.png", plot = last_plot())  # Save plot
```

---

## ✅ 7. HELPFUL `AMR` UTILITIES
```r
count_R(df$CIP_rsi)     # Number of Resistant results
count_SI(df$CIP_rsi)    # Count Susceptible + Intermediate
mo_genus("E. coli")     # Genus name
ab_atc("CIP")            # ATC code
```

---

## 🔬 8. ADVANCED: CLSI/EUCAST RULES
```r
interpretive_reading(df)         # Apply intrinsic resistance rules
mo_uncertainties(df$species)     # Detect ambiguous names
```

## 🧪 9. MDR Statistics
```r
df$mdro_result <- mdro(df, guideline = "EUCAST")
table(df$mdro_result)
```

---

## 🔍 10. LEARNING & RESOURCES
- `?AMR` → Help page for full reference
- [AMR Package Website](https://msberends.github.io/AMR/)
- [EUCAST Guidelines](https://www.eucast.org/clinical_breakpoints/)
- Kaggle: Search "AMR" under Datasets/Code
