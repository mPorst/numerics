#include <vector>

std::vector<double> euler( double h, double N0, double tau1, double tau2 );
std::vector<double> analytical( double h, double N0, double tau1, double tau2 );
std::vector<double> eulerModified( double h, double N0, double tau1, double tau2 );
std::vector<double> RK3a( double h, double N0, double tau1, double tau2 );
std::vector<double> verlet( double h, double N0, double tau1, double tau2 );
