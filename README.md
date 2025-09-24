# Cox-Ross-Rubinstein Binomial Option Pricing Model And Leisen-Reimer Implementation

This is a C++ implementation of the Cox-Ross-Rubinstein (CRR) binomial option pricing model. The model is used for pricing options by constructing a binomial lattice representing possible paths of the underlying asset price.

## Features

- Pricing of European and American options (both calls and puts)
- Calculation of option Greeks (Delta, Gamma, Theta, Vega, Rho)
- Flexible number of time steps in the binomial tree
- Interactive command-line interface

## Mathematical Model

The CRR model uses the following parameters:
- `u`: Up factor (price increase when moving up in the tree)
- `d`: Down factor (price decrease when moving down in the tree)
- `p`: Risk-neutral probability of an up move

These parameters are calculated as:
- `u = exp(σ * sqrt(Δt))`
- `d = 1/u`
- `p = (exp(r * Δt) - d) / (u - d)`

Where:
- `σ` is the volatility
- `r` is the risk-free interest rate
- `Δt` is the time step size (T/steps)

## Building the Code

Compile the program using a C++ compiler:

```bash
g++ -o crr cox_ross_rubinstein.cpp -std=c++11
```

## Running the Program

Execute the compiled program:

```bash
./crr
```

You will be prompted to enter the following parameters:
- Initial stock price
- Strike price
- Risk-free rate (as a decimal, e.g., 0.05 for 5%)
- Volatility (as a decimal, e.g., 0.2 for 20%)
- Time to maturity (in years)
- Number of steps in the binomial tree
- Option type (call/put)
- Exercise type (european/american)

## Example Usage

```
Enter initial stock price: 100
Enter strike price: 100
Enter risk-free rate (decimal, e.g., 0.05 for 5%): 0.05
Enter volatility (decimal, e.g., 0.2 for 20%): 0.2
Enter time to maturity (in years): 1
Enter number of steps in the binomial tree: 100
Enter option type (call/put): call
Enter exercise type (european/american): european

Cox-Ross-Rubinstein Model Parameters:
Up factor (u): 1.02020
Down factor (d): 0.98020
Risk-neutral probability (p): 0.51208

Option Price: 10.450383
Delta: 0.624756
Gamma: 0.017246
Theta (per day): -0.012005
Vega (for 1% change): 0.392710
Rho (for 1% change): 0.486953
```

## Implementation Details

The implementation uses a triangular array to store the stock and option price trees. The index of a node at step `i` and position `j` is calculated as `(i * (i + 1)) / 2 + j`.

### Tree Structure

For a 3-step tree:
```
       [0]
     /    \
   [1]     [2]
  /  \    /  \
[3]   [4] [5]  [6]
```

The indexes are mapped as follows:
- Step 0: [0]
- Step 1: [1], [2]
- Step 2: [3], [4], [5] 
- Step 3: [6], [7], [8], [9]

### Option Greeks

The Greeks are calculated using finite difference methods:
- Delta: First derivative of option price with respect to underlying price
- Gamma: Second derivative of option price with respect to underlying price
- Theta: Rate of change of option price with respect to time
- Vega: Rate of change of option price with respect to volatility
- Rho: Rate of change of option price with respect to interest rate 

## Leisen-Reimer Binomial Tree (lrtree.cpp)

This file implements the Leisen-Reimer binomial lattice model (1995), which provides smooth and fast convergence to Black-Scholes prices compared to the standard Cox-Ross-Rubinstein model.

### Features

- Pricing of European and American options (both calls and puts)
- Superior convergence properties compared to CRR model
- Calculation of option Greeks using finite difference methods
- Automatic adjustment to odd number of steps (required by LR method)
- Interactive command-line interface

### Mathematical Model

The Leisen-Reimer model uses the Peizer-Pratt inversion of the normal cumulative distribution function to determine the lattice parameters:

#### Key Parameters:
- **u**: Up factor = `drift * (p' / p)`
- **d**: Down factor = `drift * ((1 - p') / (1 - p))`
- **p**: Risk-neutral probability = `PP_inv(d2, N)`
- **p'**: Modified probability = `PP_inv(d1, N)`

Where:
- `drift = exp(r * Δt)`
- `d1 = (ln(S/K) + (r + 0.5*σ²)*T) / (σ*sqrt(T))`
- `d2 = d1 - σ*sqrt(T)`
- `PP_inv(z, n)` is the Peizer-Pratt inversion function

#### Peizer-Pratt Inversion Formula:
```
PP_inv(z, n) = 0.5 + 0.5 * sign(z) * sqrt(1 - exp(-z²*b/a²))
```
Where:
- `a = n + 1/3 + 0.1/(n+1)`
- `b = n + 1/6`

### Building

```bash
g++ -o lrtree lrtree.cpp -std=c++11
```

### Running

```bash
./lrtree
```

You will be prompted to enter:
- Spot price
- Strike price
- Risk-free rate (as decimal, e.g., 0.05 for 5%)
- Volatility (as decimal, e.g., 0.2 for 20%)
- Time to maturity (in years)
- Number of steps (odd numbers preferred, automatically adjusted if even)
- Option type (call/put)
- Exercise type (european/american)

### Implementation Details

The LRTree class uses backward induction through the binomial lattice:

1. **Calibration**: Automatically ensures odd number of steps and calculates u, d, p using Leisen-Reimer formulas
2. **Pricing**: Uses standard backward induction with early exercise checks for American options
3. **Greeks**: Calculated using bump-and-reprice finite difference methods

### Advantages over CRR Model

- **Faster Convergence**: Requires fewer steps to achieve accurate Black-Scholes prices
- **Smooth Convergence**: Eliminates oscillations common in CRR model
- **Better Accuracy**: Particularly effective for at-the-money options

## Black-Scholes Model (black_scholes.cpp)

This file implements the Black-Scholes formula for pricing European options and calculating their associated Greeks.

### Features

- Pricing of European call and put options.
- Calculation of Greeks: Delta, Gamma, Theta, Vega, Rho.
- Interactive command-line interface.

### Mathematical Formulas

- **d1** = `(ln(S/K) + (r + 0.5 * σ^2) * T) / (σ * sqrt(T))`
- **d2** = `d1 - σ * sqrt(T)`
- **Call Price** = `S * N(d1) - K * e^(-rT) * N(d2)`
- **Put Price** = `K * e^(-rT) * N(-d2) - S * N(-d1)`

Where `N(x)` is the standard normal cumulative distribution function.

### Building

```bash
g++ -o bs black_scholes.cpp -std=c++11 -lm
```
*(Note: `-lm` might be needed on some Linux systems to link the math library)*

### Running

```bash
./bs
```
The program will prompt for the necessary input parameters (S, K, T, r, σ).

### Greeks Calculation Details

The Greeks are calculated using their standard Black-Scholes formulas.
- **Theta** is reported per day (annualized theta / 365).
- **Vega** is reported per 1% change in volatility (calculated vega / 100).
- **Rho** is reported per 1% change in the risk-free rate (calculated rho / 100).

## Contact
[lmrojas99@gmail.com]
https://www.linkedin.com/in/lmrojasv
