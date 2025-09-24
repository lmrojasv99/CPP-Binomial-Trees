#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

/*
   Leisen–Reimer binomial lattice (1995)
   --------------------------------------------------
   * Smooth and fast convergence to Black–Scholes.
   * Requires an odd number of steps (internally rounded up).
   * Up/Down factors and risk–neutral probability follow the
     Peizer–Pratt inversion of the normal CDF as described in
     Leisen & Reimer (1995) and many secondary sources.
*/

class LRTree {
private:
    // Inputs
    double S0;      // spot
    double K;       // strike
    double r;       // risk-free rate (continuously compounded)
    double sigma;   // volatility
    double T;       // maturity (years)
    int    N;       // steps (forced to odd)
    bool   isCall;  // true = call, false = put
    bool   isAmerican; // true = American, false = European

    // Lattice parameters
    double dt{},    // time step
           u{},     // up factor
           d{},     // down factor
           p{};     // risk-neutral probability

    /* Black–Scholes d1 / d2 */
    static double d1(double S, double K, double r, double sigma, double tau) {
        return (std::log(S / K) + (r + 0.5 * sigma * sigma) * tau) /
               (sigma * std::sqrt(tau));
    }
    static double d2(double S, double K, double r, double sigma, double tau) {
        return d1(S, K, r, sigma, tau) - sigma * std::sqrt(tau);
    }

    /* Peizer–Pratt inversion   h^{-1}(z)  */
    static double ppInversion(double z, int n) {
        if (z == 0.0) return 0.5;           // symmetry
        double a = n + 1.0 / 3.0 + 0.1 / (n + 1.0);
        double b = n + 1.0 / 6.0;
        double x = z / a;
        double expTerm = std::exp(-x * x * b);
        double sgn = z > 0.0 ? 1.0 : -1.0;
        return 0.5 + 0.5 * sgn * std::sqrt(1.0 - expTerm);
    }

    /* Derive u, d, p according to Leisen–Reimer */
    void calibrate() {
        // ensure odd number of steps as required by LR
        if (N % 2 == 0) ++N;
        dt = T / static_cast<double>(N);

        double d_1 = d1(S0, K, r, sigma, T);
        double d_2 = d_1 - sigma * std::sqrt(T);

        double pPrime = ppInversion(d_1, N);
        p             = ppInversion(d_2, N);

        // cost of carry b = r (no dividend); extend easily if q given
        double drift = std::exp(r * dt);

        u = drift * (pPrime / p);
        d = drift * ((1.0 - pPrime) / (1.0 - p));
    }

    /* helper: payoff */
    inline double payoff(double ST) const {
        if (isCall) return std::max(0.0, ST - K);
        return std::max(0.0, K - ST);
    }

public:
    LRTree(double S0_, double K_, double r_, double sigma_, double T_, int steps,
           bool callOption, bool american)
        : S0(S0_), K(K_), r(r_), sigma(sigma_), T(T_), N(steps),
          isCall(callOption), isAmerican(american) {
        calibrate();
    }

    // Price via backward induction
    double price() const {
        std::vector<double> option(N + 1);

        // terminal values
        for (int i = 0; i <= N; ++i) {
            double ST = S0 * std::pow(u, i) * std::pow(d, N - i);
            option[i] = payoff(ST);
        }

        double disc = std::exp(-r * dt);

        // backward pass
        for (int step = N - 1; step >= 0; --step) {
            for (int i = 0; i <= step; ++i) {
                double continuation = disc * (p * option[i + 1] + (1.0 - p) * option[i]);
                if (isAmerican) {
                    double ST = S0 * std::pow(u, i) * std::pow(d, step - i);
                    option[i] = std::max(continuation, payoff(ST));
                } else {
                    option[i] = continuation;
                }
            }
        }
        return option[0];
    }

    // Bump-and-reprice Greeks (coarse but simple)
    double delta(double dS = 0.01) {
        LRTree up(S0 + dS, K, r, sigma, T, N, isCall, isAmerican);
        LRTree down(S0 - dS, K, r, sigma, T, N, isCall, isAmerican);
        return (up.price() - down.price()) / (2.0 * dS);
    }
    double gamma(double dS = 0.01) {
        LRTree up(S0 + dS, K, r, sigma, T, N, isCall, isAmerican);
        LRTree mid(S0,       K, r, sigma, T, N, isCall, isAmerican);
        LRTree down(S0 - dS, K, r, sigma, T, N, isCall, isAmerican);
        return (up.price() - 2.0 * mid.price() + down.price()) / (dS * dS);
    }
    double vega(double dVol = 0.01) {
        LRTree highVol(S0, K, r, sigma + dVol, T, N, isCall, isAmerican);
        LRTree lowVol (S0, K, r, sigma - dVol, T, N, isCall, isAmerican);
        return (highVol.price() - lowVol.price()) / (2.0 * dVol);
    }
    double theta(double dT = 1.0 / 365.0) {
        if (T <= dT) return 0.0;
        LRTree lessT(S0, K, r, sigma, T - dT, N, isCall, isAmerican);
        // Market convention: theta = dV/dt = -dV/dτ where τ is time-to-maturity.
        // Since (lessT.price() - price()) ≈ -dV, divide by +dT to obtain the correct sign.
        return (lessT.price() - price()) / dT; // per year, calendar-time decay
    }
    double rho(double dR = 0.0001) {
        LRTree highR(S0, K, r + dR, sigma, T, N, isCall, isAmerican);
        LRTree lowR (S0, K, r - dR, sigma, T, N, isCall, isAmerican);
        return (highR.price() - lowR.price()) / (2.0 * dR);
    }
};

int main() {
    std::cout << "Leisen–Reimer Binomial Tree Option Pricer\n";

    double S0, K, r, sigma, T;
    int steps;
    std::string optType, exercise;

    std::cout << "Spot price: ";                       std::cin >> S0;
    std::cout << "Strike price: ";                     std::cin >> K;
    std::cout << "Risk-free rate (e.g. 0.05): ";       std::cin >> r;
    std::cout << "Volatility (e.g. 0.2): ";            std::cin >> sigma;
    std::cout << "Time to maturity in years: ";        std::cin >> T;
    std::cout << "Number of steps (odd preferred): ";  std::cin >> steps;
    std::cout << "Option type (call/put): ";           std::cin >> optType;
    std::cout << "Exercise (european/american): ";     std::cin >> exercise;

    
    std::transform(optType.begin(), optType.end(), optType.begin(), ::tolower);
    std::transform(exercise.begin(), exercise.end(), exercise.begin(), ::tolower);

    bool isCall   = (optType == "call" || optType == "c");
    bool american = (exercise == "american" || exercise == "a");

    LRTree pricer(S0, K, r, sigma, T, steps, isCall, american);

    double optionPrice = pricer.price();
    double Delta = pricer.delta();
    double Gamma = pricer.gamma();
    double Vega  = pricer.vega();
    double Theta = pricer.theta();
    double Rho   = pricer.rho();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nOption price: " << optionPrice << "\n";
    std::cout << "Delta:  " << Delta  << "\n";
    std::cout << "Gamma:  " << Gamma  << "\n";
    std::cout << "Vega:   " << Vega   << "\n";
    std::cout << "Theta:  " << Theta  << "\n";
    std::cout << "Rho:    " << Rho    << "\n";

    return 0;
}
