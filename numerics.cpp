#include "numerics.h"
#include <vector>
#include <math.h>
#include <iostream>

// equation
double equation::f(double t, unsigned int iterStep,std::vector<double> y)
{
	// t is unused because the diff. eq. for radioactive decay does not contain t explicitly. Not the case in general !
	return -l1*y.at(iterStep-1); 
}
		
double equation::gety0() {return y0;}


// solver
std::vector<double> solver::gety() {return y;}
std::vector<double> solver::gett() {return t;}
unsigned int solver::getIterStep() {return iterStep;}

void solver::resetSolver()
{
	iterStep=0;
	y.erase(y.begin(), y.end());
	t.erase(t.begin(), t.end());
}

// euler

euler::euler() 
{iterStep = 0;}


double euler::iterateY(double h) 
{
	t.push_back(iterStep*h);

	if(iterStep > 0){
		y.push_back( y.at(iterStep-1)+h*solveEq.f(t.at(iterStep-1), iterStep, y) ); // t.at(iterStep) or t.at(iterStep-1) ???
	}
	else{
		y.push_back( solveEq.gety0());
	}
	return y.at(iterStep++); // the whole y vector gets stored internally - so return of the iterated y is a bonus
}

void euler::iterateYN(double h, unsigned int n)
{
	for(unsigned int i=0; i<n; i++){
		this->iterateY(h);
	}
}


std::vector<double> analytical( double h, double N0, double tau1, double tau2 )
{
	double t;
	std::vector<double> N;
	for(int i=0; i<1000; i++)
	{
		t=i*h;
		N.push_back(N0*exp(-1./tau1*t));
	}
	return N;
}

std::vector<double> eulerf( double h, double N0, double tau1, double tau2 )
{
	double t;
	std::vector<double> N;
	for(int i=0; i<1000; i++)
	{
		t=i*h;

		if(i>0){
		N.push_back( N[i-1]-1./(tau1)*N[i-1]*h );
		//std::cout << N[i] << std::endl;
		}
		else{
		N.push_back(N0);
		}
	}
	return N;
}

std::vector<double> eulerModified( double h, double N0, double tau1, double tau2 )
{
double t;
std::vector<double> N;
for(int i=0; i<1000; i++)
{
	t=i*h;

	if(i>0){
	N.push_back(N[i-1]-(h*(1./tau1*N[i-1]+1./tau1*(N[i-1]-1./tau1*h*N[i-1]))/2) );
	//std::cout << N[i-1]*h*1./tau1 << '\t' << h*pow(1./tau1,2)*N[i-1] << std::endl;
	}
	else{
	N.push_back(N0);
	}
}
return N;
}

std::vector<double> RK3a( double h, double N0, double tau1, double tau2 )
{
std::vector<double> N;
double l1 = 1./tau1;

double rkpart1, rkpart2, rkpart3;

for(int i=0; i<1000; i++)
{
	//double t=i*h;
	if(i>0){
	rkpart1 = -h*l1*N.at(i-1);
	rkpart2 = -4*h*l1*(N.at(i-1)-h*l1*N.at(i-1)/2);
	rkpart3 = -h*l1*(N.at(i-1)+h*l1*N.at(i-1)-2*h*l1*(N.at(i-1)-h*l1*N.at(i-1)/2));
	N.push_back(N.at(i-1)+(1./6)*(rkpart1+rkpart2+rkpart3)); 
//	std::cout << N.at(i-1) << rkpart1 << '\t' << rkpart2 << '\t' << rkpart3 << '\t' << std::endl;
	}
	else {
	N.push_back(N0);
	}
}
return N;
}	

std::vector<double> verlet( double h, double N0, double tau1, double tau2 )
{
	std::vector<double> N;
	double l1 = 1./tau1;

	for(int i=0; i<1000; i++)
	{
		if(i>1){
			N.push_back(2*N[i-1]-N[i-2]+h*h*l1*l1*N[i-1]);
			//std::cout << "i: " << i << '\t' << "h*h*l1*l1*N[i-1]: " << h*h*l1*l1*N[i-1] << std::endl;
			//std::cout << "i: " << i << '\t' << "2*Ni-1: " << 2*N[i-1] << std::endl;
			//std::cout << "i: " << i << '\t' << "Ni-2: " << N[i-2] << std::endl;
		}
		else if(i==1) {
			N.push_back(N0-h*l1*N0);
		}
		else{ // if i==0
			N.push_back(N0);
		}
	}
	return N;
			
}
