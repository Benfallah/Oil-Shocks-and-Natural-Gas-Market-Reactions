# Oil Shocks and Natural Gas Market Reactions

This project is my replication and exploration of how the U.S. natural gas market responds to demand and supply shocks in the global crude oil market.  
It uses **VAR models** to trace shocks, compute impulse responses, and look at cumulative contributions of different oil shocks on natural gas prices.

---

## What the Code Does
- **Data preparation**: Reads in monthly time series for oil production, real activity, crude oil price, and natural gas price.  
- **VAR estimation**: Fits a vector autoregression to capture the dynamics between these variables.  
- **Structural shocks**: Uses a Cholesky decomposition to identify supply, demand, and oil price shocks.  
- **Impulse responses (IRFs)**: Traces the effect of each shock over time.  
- **Cumulative effects**: Adds up how much each shock explains movements in oil prices.  
- **Plots and tables**: Outputs time series charts, impulse response figures, and summary tables for easier interpretation.

---

## Why This Matters
Understanding how **oil market shocks spill over into natural gas** is important for:
- Policy analysis (energy security, climate policy)
- Market participants (hedging and forecasting)
- Researchers interested in energy economics

This is also a personal learning project, showing my workflow in **R (tidyverse + vars package)** with reproducible coding practices.

---

## Repo Structure
```
.
├── Codes_Oil_Shocks.R          # main script (original)
├── Oil_Shocks_VAR_Refactor.R   # cleaned and refactored version (tidyverse style)
├── data/                       # input data (not uploaded here)
├── output/                     # optional folder for plots/tables
└── README.md                   # this file
```

---

## How to Run
1. Clone this repo.  
2. Place the dataset (`data.txt`) inside a `data/` folder.  
   - Expected columns: `prod`, `ea`, `po`, `pg`  
   - Monthly frequency, starting in 1978.  
3. Run `Oil_Shocks_VAR_Refactor.R` in RStudio (R 4.0+ recommended).  
4. Outputs: plots of series, IRFs, and cumulative contributions.

---


---

## Personal Note
This repo is part of my **energy economics portfolio**, where I explore links between oil and gas markets using econometric tools.  
It reflects both my research interests and my effort to write clean, reproducible R code.
