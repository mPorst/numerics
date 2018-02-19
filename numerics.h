#include <vector>

class equation{
	private:
		double tau1 = 100000;
		double l1 = 1./tau1;
		double y0=100000000;

	public:
		double f(double t, double y); // for radioactive decay: -l1*y[i-1]
		double gety0();
};

class solver{
	protected:
		std::vector<double> y; // x'(t) = f(t,x(t)) and y(t) are the numerical approximations
		std::vector<double> t;
		equation solveEq;
		unsigned int iterStep=0;

	public:
		virtual void iterateY(double h) = 0;
		void iterateYN(double h, unsigned int n);
		std::vector<double> gety();
		std::vector<double> gett();
		unsigned int getIterStep();
		void resetSolver();
};

class euler: public solver{
	public:
		void iterateY(double h);		
};

class eulerModified: public solver{
	public:
		void iterateY(double h);		
};

class RK3a: public solver{ // runge kutta 3rd order variant
	public:
		void iterateY(double h);		
};


std::vector<double> analytical( double h, double N0, double tau1, double tau2 );
std::vector<double> verlet( double h, double N0, double tau1, double tau2 );
