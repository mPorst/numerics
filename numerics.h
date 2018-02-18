#include <vector>

class equation{
	private:
		double tau1 = 100000;
		double l1 = 1./tau1;
		double y0=100000000;

	public:
		double f(double t, unsigned int iterStep,std::vector<double> y); // for radioactive decay: -l1*y[i-1]
		double gety0();
};

class solver{
	protected:
		std::vector<double> y; // x'(t) = f(t,x(t)) and y(t) are the numerical approximations
		std::vector<double> t;
		equation solveEq;
		unsigned int iterStep;

	public:
		virtual double iterateY(double h) = 0;
		virtual void iterateYN(double h, unsigned int n) = 0;
		std::vector<double> gety();
		std::vector<double> gett();
		unsigned int getIterStep();
		void resetSolver();
};

class euler: public solver{
	public:
		euler();
		double iterateY(double h);		
		void iterateYN(double h, unsigned int n);
};



std::vector<double> eulerf( double h, double N0, double tau1, double tau2 );
std::vector<double> analytical( double h, double N0, double tau1, double tau2 );
std::vector<double> eulerModified( double h, double N0, double tau1, double tau2 );
std::vector<double> RK3a( double h, double N0, double tau1, double tau2 );
std::vector<double> verlet( double h, double N0, double tau1, double tau2 );
