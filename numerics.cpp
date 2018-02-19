#include "numerics.h"
#include <vector>
#include <math.h>
#include <iostream>

// equation
double equation::f(double t, double y)
{
	// t is unused because the diff. eq. for radioactive decay does not contain t explicitly. Not the case in general !
	return -l1*y; 
}
		
double equation::gety0() {return y0;}


// solver
std::vector<double> solver::gety() {return y;}
std::vector<double> solver::gett() {return t;}
unsigned int solver::getIterStep() {return iterStep;}

void solver::iterateYN(double h, unsigned int n)
{
	for(unsigned int i=0; i<n; i++){
		this->iterateY( h );
	}
}

void solver::resetSolver()
{
	iterStep=0;
	y.erase(y.begin(), y.end());
	t.erase(t.begin(), t.end());
}

// euler

//euler::euler() 
//{iterStep = 0;}


void euler::iterateY(double h) 
{
	t.push_back(iterStep*h);

	if(iterStep > 0){
		y.push_back( y.at(iterStep-1)+h*solveEq.f(t.at(iterStep-1), y.at(iterStep-1)) ); // t.at(iterStep) or t.at(iterStep-1) ???
	}
	else{
		y.push_back( solveEq.gety0());
	}
	++iterStep;
}


// euler modified

void eulerModified::iterateY(double h) 
{
	t.push_back(iterStep*h);

	if(iterStep > 0){
		y.push_back(y.at(iterStep-1)+0.5*h*(solveEq.f(t.at(iterStep-1), y.at(iterStep-1))+solveEq.f(t.at(iterStep), y.at(iterStep-1)+h*solveEq.f(t.at(iterStep-1), y.at(iterStep-1)))));
	}
	else{
		y.push_back( solveEq.gety0());
	}
	++iterStep;
}

// rk3a

void RK3a::iterateY(double h)
{
	t.push_back(iterStep*h);

	if(iterStep > 0){
		y.push_back(y.at(iterStep-1)+1./6*(h*solveEq.f(t.at(iterStep-1), y.at(iterStep-1))+4*h*solveEq.f(t.at(iterStep-1)+h/2, y.at(iterStep-1)+h*solveEq.f(t.at(iterStep-1),y.at(iterStep-1))/2)+h*solveEq.f(t.at(iterStep-1)+h, y.at(iterStep-1)-h*solveEq.f(t.at(iterStep-1), y.at(iterStep-1))+2*h*solveEq.f(t.at(iterStep-1)-h/2, y.at(iterStep-1)+h*solveEq.f(t.at(iterStep-1), y.at(iterStep-1))/2))));
	}
	else{
		y.push_back( solveEq.gety0());
	}
	++iterStep;
}


std::vector<double> analytical( double h, double N0, double tau1, double tau2 )
{
	double t;
	std::vector<double> N;
	for(int i=0; i<1000; i++){
		t=i*h;
		N.push_back(N0*exp(-1./tau1*t));
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
