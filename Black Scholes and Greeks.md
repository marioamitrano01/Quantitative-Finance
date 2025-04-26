# Black-Scholes Option Pricing Model in Julia

A comprehensive Julia implementation of the Black-Scholes model for pricing European options and calculating Greeks.

## Features

- Price European call and put options
- Calculate all option Greeks (Delta, Gamma, Vega, Theta, Rho)
- Generate high-quality visualizations
- Verify mathematical relationships (put-call parity)
- Individual plots for each Greek

## Requirements

```julia
import Pkg
Pkg.add.(["Plots", "Distributions", "Dates"])
```
## Output

The analysis creates a folder on your Desktop containing:
- Price surface plots for both call and put options
- Individual plots for each Greek (Delta, Gamma, Vega, Theta, Rho)
- Price vs. volatility and time decay plots
- Summary text files with numerical results and formulas
