#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

// Define M_PI if not defined (e.g., on Windows with MSVC)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class BlackScholesModel {
private:
    double S;
    double K;
    double T;
    double r;
    double sigma;

    // Helper function for Standard Normal CDF using erf
    double norm_cdf(double value) const {
        return 0.5 * std::erfc(-value * M_SQRT1_2); // M_SQRT1_2 is 1/sqrt(2)
    }

    // Helper function for Standard Normal PDF
    double norm_pdf(double value) const {
        return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * value * value);
    }

    // Calculate d1
    double calculate_d1() const {
        if (sigma == 0 || T == 0) return std::numeric_limits<double>::infinity(); // Avoid division by zero
        return (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    }

    // Calculate d2
    double calculate_d2(double d1_val) const {
         if (sigma == 0 || T == 0) return std::numeric_limits<double>::infinity(); // Avoid division by zero
        return d1_val - sigma * std::sqrt(T);
    }

public:
    // Constructor
    BlackScholesModel(double stockPrice, double strikePrice, double timeToExp, double riskFreeRate, double volatility)
        : S(stockPrice), K(strikePrice), T(timeToExp), r(riskFreeRate), sigma(volatility) {}

    double call_option_price() const {
        double d1_val = calculate_d1();
        // Handle edge cases or invalid inputs for d1
        if (!std::isfinite(d1_val)) return (S > K * exp(-r*T)) ? S - K * exp(-r*T) : 0.0; // Intrinsic value or 0
        double d2_val = calculate_d2(d1_val);
        return S * norm_cdf(d1_val) - K * std::exp(-r * T) * norm_cdf(d2_val);
    }

    double put_option_price() const {
        double d1_val = calculate_d1();
         // Handle edge cases or invalid inputs for d1
        if (!std::isfinite(d1_val)) return (K * exp(-r*T) > S) ? K * exp(-r*T) - S : 0.0; // Intrinsic value or 0
        double d2_val = calculate_d2(d1_val);
        return K * std::exp(-r * T) * norm_cdf(-d2_val) - S * norm_cdf(-d1_val);
    }

    double delta_call() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val)) return (S >= K * exp(-r*T)) ? 1.0 : 0.0;
        return norm_cdf(d1_val);
    }

    double delta_put() const {
        double d1_val = calculate_d1();
         if (!std::isfinite(d1_val)) return (S < K * exp(-r*T)) ? -1.0 : 0.0;
        return norm_cdf(d1_val) - 1.0; // delta_put = delta_call - 1
    }

    double gamma() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val) || S <= 0 || sigma <= 0 || T <= 0) return 0.0;
        return norm_pdf(d1_val) / (S * sigma * std::sqrt(T));
    }

    double theta_call() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val) || S <= 0 || sigma <= 0 || T <= 0) return 0.0;
        double d2_val = calculate_d2(d1_val);
        double term1 = - (S * norm_pdf(d1_val) * sigma) / (2.0 * std::sqrt(T));
        double term2 = - r * K * std::exp(-r * T) * norm_cdf(d2_val);
        return term1 + term2;
    }

     double theta_put() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val) || S <= 0 || sigma <= 0 || T <= 0) return 0.0;
        double d2_val = calculate_d2(d1_val);
        double term1 = - (S * norm_pdf(d1_val) * sigma) / (2.0 * std::sqrt(T));
        double term2 = + r * K * std::exp(-r * T) * norm_cdf(-d2_val); // Note the + sign and -d2_val
        return term1 + term2;
    }

    double vega() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val) || S <= 0 || T <= 0) return 0.0;
        return S * norm_pdf(d1_val) * std::sqrt(T);
    }

    double rho_call() const {
        double d1_val = calculate_d1();
        if (!std::isfinite(d1_val) || T <= 0) return 0.0;
        double d2_val = calculate_d2(d1_val);
        return K * T * std::exp(-r * T) * norm_cdf(d2_val);
    }

    double rho_put() const {
        double d1_val = calculate_d1();
         if (!std::isfinite(d1_val) || T <= 0) return 0.0;
        double d2_val = calculate_d2(d1_val);
        return -K * T * std::exp(-r * T) * norm_cdf(-d2_val);
    }
};

int main() {
    // Variables for user input
    double S, K, T, r, sigma;

    // Get user input
    std::cout << "--- Black-Scholes Model ---" << std::endl;
    std::cout << "Enter underlying asset price (S): ";
    std::cin >> S;
    std::cout << "Enter option strike price (K): ";
    std::cin >> K;
    std::cout << "Enter time to expiration in years (T): ";
    std::cin >> T;
    std::cout << "Enter risk-free interest rate (r, e.g., 0.05 for 5%): ";
    std::cin >> r;
    std::cout << "Enter volatility (sigma, e.g., 0.2 for 20%): ";
    std::cin >> sigma;

    // Validate inputs
    if (S <= 0 || K <= 0 || T <= 0 || sigma < 0) {
         std::cerr << "Error: Invalid input parameters. Prices, time, and strike must be positive. Volatility cannot be negative." << std::endl;
         return 1;
    }
     if (sigma == 0) {
        std::cout << "\nWarning: Volatility is zero. Greeks calculations involving sigma may yield infinite or NaN results or simplified intrinsic value." << std::endl;
    }


    // Create Black-Scholes model instance
    BlackScholesModel bsm(S, K, T, r, sigma);

    // Calculate and display results
    std::cout << std::fixed << std::setprecision(6); // Set output precision

    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "Call Option Price: " << bsm.call_option_price() << std::endl;
    std::cout << "Put Option Price:  " << bsm.put_option_price() << std::endl;

    std::cout << "\n--- Greeks ---" << std::endl;
    std::cout << "Delta (Call):      " << bsm.delta_call() << std::endl;
    std::cout << "Delta (Put):       " << bsm.delta_put() << std::endl;
    std::cout << "Gamma:             " << bsm.gamma() << std::endl;
    // Theta is typically quoted per trading day (252 trading days per year in US markets)
    std::cout << "Theta (Call, /day):" << bsm.theta_call() / 252.0 << std::endl;
    std::cout << "Theta (Put, /day): " << bsm.theta_put() / 252.0 << std::endl;
    // Vega is typically quoted per 1% change in volatility (0.01)
    std::cout << "Vega (per 1% vol): " << bsm.vega() / 100.0 << std::endl;
     // Rho is typically quoted per 1% change in interest rate (0.01)
    std::cout << "Rho (Call, per 1% rate): " << bsm.rho_call() / 100.0 << std::endl;
    std::cout << "Rho (Put, per 1% rate):  " << bsm.rho_put() / 100.0 << std::endl;

    return 0;
} 