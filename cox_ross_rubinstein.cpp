#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>

class CoxRossRubinsteinTree {
private:
    int steps;
    double S0;          // Initial stock price
    double K;           // Strike price
    double r;           // Risk-free rate
    double sigma;       // Volatility
    double T;           // Time to maturity
    bool isCall;        // Option type (call or put)
    bool isAmerican;    // Exercise type (American or European)
    
    std::vector<double> stockTree;   // Stock price tree
    std::vector<double> optionTree;  // Option price tree
    
    // CRR parameters
    double u;           // Up factor
    double d;           // Down factor
    double p;           // Risk-neutral probability
    
    // Calculate CRR parameters
    void calculateParameters() {
        double dt = T / steps;
        
        // Cox-Ross-Rubinstein parameters
        u = exp(sigma * sqrt(dt));
        d = 1.0 / u;
        p = (exp(r * dt) - d) / (u - d);
    }

public:
    CoxRossRubinsteinTree(double initialPrice, double strikePrice, double riskFreeRate, 
                   double volatility, double timeToMaturity, int numberOfSteps, 
                   bool callOption = true, bool americanOption = false) 
        : S0(initialPrice), K(strikePrice), r(riskFreeRate), sigma(volatility),
          T(timeToMaturity), steps(numberOfSteps), isCall(callOption), isAmerican(americanOption) {
        
        // Calculate the model parameters
        calculateParameters();
        
        // Initialize trees
        int treeSize = (steps + 1) * (steps + 2) / 2; // Triangular number
        stockTree.resize(treeSize);
        optionTree.resize(treeSize);
    }
    
    // Get node index in the tree
    int getIndex(int step, int node) const {
        return (step * (step + 1)) / 2 + node;
    }
    
    // Calculate stock price at a specific node
    double getStockPrice(int step, int node) const {
        return S0 * pow(u, node) * pow(d, step - node);
    }
    
    // Calculate option payoff at expiration
    double calculatePayoff(double stockPrice) const {
        if (isCall) {
            return std::max(0.0, stockPrice - K);
        } else {
            return std::max(0.0, K - stockPrice);
        }
    }
    
    // Calculate option price
    double price() {
        // Build the stock price tree
        for (int i = 0; i <= steps; i++) {
            for (int j = 0; j <= i; j++) {
                int idx = getIndex(i, j);
                stockTree[idx] = getStockPrice(i, j);
            }
        }
        
        // Calculate option values at expiration (last step)
        for (int j = 0; j <= steps; j++) {
            int idx = getIndex(steps, j);
            optionTree[idx] = calculatePayoff(stockTree[idx]);
        }
        
        // Work backwards through the tree
        double dt = T / steps;
        double discountFactor = exp(-r * dt);
        
        for (int i = steps - 1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                int currentIdx = getIndex(i, j);
                int upIdx = getIndex(i + 1, j + 1);
                int downIdx = getIndex(i + 1, j);
                
                // Calculate expected value
                double expectedValue = discountFactor * (p * optionTree[upIdx] + (1 - p) * optionTree[downIdx]);
                
                if (isAmerican) {
                    // For American options, consider early exercise
                    double exerciseValue = calculatePayoff(stockTree[currentIdx]);
                    optionTree[currentIdx] = std::max(expectedValue, exerciseValue);
                } else {
                    // For European options
                    optionTree[currentIdx] = expectedValue;
                }
            }
        }
        
        // Return the option price at the root of the tree
        return optionTree[0];
    }
    
    // Calculate option Greeks
    
    // Delta - sensitivity to change in underlying price
    double calculateDelta() {
        if (steps < 1) return 0.0;
        
        // Use first level of the tree to calculate delta
        int upIdx = getIndex(1, 1);
        int downIdx = getIndex(1, 0);
        
        double upOptionValue = optionTree[upIdx];
        double downOptionValue = optionTree[downIdx];
        double upStockPrice = stockTree[upIdx];
        double downStockPrice = stockTree[downIdx];
        
        return (upOptionValue - downOptionValue) / (upStockPrice - downStockPrice);
    }
    
    // Gamma - sensitivity of delta to change in underlying price
    double calculateGamma() {
        if (steps < 2) return 0.0;
        
        // First, price the option to build the tree
        if (optionTree[0] == 0) price();
        
        // Calculate delta at up and down nodes at step 1
        double deltaUp, deltaDown;
        
        // Up node delta
        int upUpIdx = getIndex(2, 2);
        int upDownIdx = getIndex(2, 1);
        deltaUp = (optionTree[upUpIdx] - optionTree[upDownIdx]) / 
                 (stockTree[upUpIdx] - stockTree[upDownIdx]);
        
        // Down node delta
        int downUpIdx = getIndex(2, 1); // Same as upDownIdx
        int downDownIdx = getIndex(2, 0);
        deltaDown = (optionTree[downUpIdx] - optionTree[downDownIdx]) / 
                   (stockTree[downUpIdx] - stockTree[downDownIdx]);
        
        // Stock prices at step 1
        int upIdx = getIndex(1, 1);
        int downIdx = getIndex(1, 0);
        double upStockPrice = stockTree[upIdx];
        double downStockPrice = stockTree[downIdx];
        
        // Calculate gamma
        return (deltaUp - deltaDown) / (upStockPrice - downStockPrice);
    }
    
    // Theta - sensitivity to change in time to expiration
    double calculateTheta() {
        // Approximate theta by creating a new tree with slightly less time
        double originalPrice = price();
        
        // Create a new tree with one day less to expiration
        double dt = 1.0 / 365.0;
        if (T <= dt) return 0.0; // Cannot reduce time any further
        
        CoxRossRubinsteinTree less_time(S0, K, r, sigma, T - dt, steps, isCall, isAmerican);
        double newPrice = less_time.price();
        
        // Calculate theta (per day)
        return (newPrice - originalPrice) / dt;
    }
    
    // Vega - sensitivity to change in volatility
    double calculateVega() {
        double originalPrice = price();
        
        // Create a new tree with slightly higher volatility
        double dv = 0.01; // 1% change in volatility
        CoxRossRubinsteinTree higher_vol(S0, K, r, sigma + dv, T, steps, isCall, isAmerican);
        double newPrice = higher_vol.price();
        
        // Calculate vega for 1% change in volatility
        return (newPrice - originalPrice) / dv;
    }
    
    // Rho - sensitivity to change in interest rate
    double calculateRho() {
        double originalPrice = price();
        
        // Create a new tree with slightly higher interest rate
        double dr = 0.01; // 1% change in rate
        CoxRossRubinsteinTree higher_rate(S0, K, r + dr, sigma, T, steps, isCall, isAmerican);
        double newPrice = higher_rate.price();
        
        // Calculate rho for 1% change in interest rate
        return (newPrice - originalPrice) / dr;
    }
    
    // Get model parameters
    double getUpFactor() const { return u; }
    double getDownFactor() const { return d; }
    double getProbability() const { return p; }
};

int main() {
    // Variables for user input
    double S0;          // Initial stock price
    double K;           // Strike price
    double r;           // Risk-free rate
    double sigma;       // Volatility
    double T;           // Time to maturity
    int steps;          // Number of steps in the tree
    std::string optionType;
    std::string exerciseType;
    
    // Get user input
    std::cout << "Enter initial stock price: ";
    std::cin >> S0;
    
    std::cout << "Enter strike price: ";
    std::cin >> K;
    
    std::cout << "Enter risk-free rate (decimal, e.g., 0.05 for 5%): ";
    std::cin >> r;
    
    std::cout << "Enter volatility (decimal, e.g., 0.2 for 20%): ";
    std::cin >> sigma;
    
    std::cout << "Enter time to maturity (in years): ";
    std::cin >> T;
    
    std::cout << "Enter number of steps in the binomial tree: ";
    std::cin >> steps;
    
    std::cout << "Enter option type (call/put): ";
    std::cin >> optionType;
    
    std::cout << "Enter exercise type (european/american): ";
    std::cin >> exerciseType;
    
    // Convert option type to boolean
    bool isCall = (optionType == "call" || optionType == "Call");
    
    // Convert exercise type to boolean
    bool isAmerican = (exerciseType == "american" || exerciseType == "American");
    
    // Create option pricing model
    CoxRossRubinsteinTree option(S0, K, r, sigma, T, steps, isCall, isAmerican);
    
    // Calculate option price and Greeks
    double price = option.price();
    double delta = option.calculateDelta();
    double gamma = option.calculateGamma();
    double theta = option.calculateTheta();
    double vega = option.calculateVega();
    double rho = option.calculateRho();
    
    // Display model parameters
    std::cout << "\nCox-Ross-Rubinstein Model Parameters:" << std::endl;
    std::cout << "Up factor (u): " << option.getUpFactor() << std::endl;
    std::cout << "Down factor (d): " << option.getDownFactor() << std::endl;
    std::cout << "Risk-neutral probability (p): " << option.getProbability() << std::endl;
    
    // Display results
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nOption Price: " << price << std::endl;
    std::cout << "Delta: " << delta << std::endl;
    std::cout << "Gamma: " << gamma << std::endl;
    std::cout << "Theta: " << theta << std::endl;
    std::cout << "Vega: " << vega << std::endl;
    std::cout << "Rho: " << rho << std::endl;
    
    return 0;
} 