#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "gnuplot-iostream.h"
#include "numerics.h"
#include <chrono>
#include "statistics.h"

//shortcut for getting time
std::chrono::high_resolution_clock::time_point getTime(){return std::chrono::high_resolution_clock::now();}

int main(){
	
	//set up gnuplot interface
	Gnuplot gp;

	//set up solver objects
	euler eulerSolver;

	double stepsize = 500;
	double numberMolecules = 100000000;
	double tau1 = 100000;
	unsigned int numberIterations = 1000;
	std::vector<double> Nana, Neuler, NeulerMod, Nrk3a, Nverlet;
	
	std::vector<std::chrono::duration<double, std::milli>> teuler, tana, teulerMod, trk3a, tverlet;
	std::fstream filestream;
	filestream.open("myfile.txt", std::fstream::out);

	// calculate 1000 times to account for variance of time measurements
	for(int i=0; i<1000; i++) {
	// call numerical methods and measure time points between them
	//
	auto t1 = getTime();
	eulerSolver.iterateYN( stepsize, numberIterations );
	Neuler = eulerSolver.gety();
	auto t2 = getTime();
	Nana = analytical(stepsize, numberMolecules, tau1, 1);
	auto t3 = getTime();
	NeulerMod = eulerModified(stepsize, numberMolecules, tau1, 1);
	auto t4 = getTime();
	Nrk3a = RK3a(stepsize, numberMolecules, tau1, 1);
	auto t5 = getTime();
	Nverlet = verlet(stepsize, numberMolecules, tau1, 1);
	auto t6 = getTime();
	
	eulerSolver.resetSolver();

	teuler.push_back(t2-t1);
	tana.push_back(t3-t2);
	teulerMod.push_back(t4-t3);
	trk3a.push_back(t5-t4);
	tverlet.push_back(t6-t5);
	}

	//mean value and standard deviation
	
	// get mean of execution time and standard deviation (defined in statistics)
	double meanteuler=0, meantana=0, meanteulerMod=0, meantrk3a=0, meantverlet=0;
	double stdteuler=0, stdtana=0, stdteulerMod=0, stdtrk3a=0, stdtverlet=0;
	
	meanteuler = mean<std::chrono::duration<double,std::milli>>(teuler);
	stdteuler = stdDev<std::chrono::duration<double,std::milli>>(teuler, meanteuler);
	meantana = mean<std::chrono::duration<double,std::milli>>(tana);
	stdtana = stdDev<std::chrono::duration<double,std::milli>>(tana, meantana);
	meanteulerMod = mean<std::chrono::duration<double,std::milli>>(teulerMod);
	stdteulerMod = stdDev<std::chrono::duration<double,std::milli>>(teulerMod, meanteulerMod);
	meantrk3a = mean<std::chrono::duration<double,std::milli>>(trk3a);
	stdtrk3a = stdDev<std::chrono::duration<double,std::milli>>(trk3a, meantrk3a);
	meantverlet = mean<std::chrono::duration<double,std::milli>>(tverlet);
	stdtverlet = stdDev<std::chrono::duration<double,std::milli>>(tverlet, meantverlet);

	std::cout << "Euler: " << meanteuler << "+-" << stdteuler << "ms"  << std::endl
		<< "analytical: " << meantana << "+-" << stdtana  << "ms" << std::endl
		<< "Runge Kutta 2nd order(mod. Euler): " << meanteulerMod << "+-" << stdteulerMod << "ms"  << std::endl
		<< "Runge Kutta 3rd order: " << meantrk3a << "+-" << stdtrk3a << "ms"  << std::endl
		<< "verlet: " << meantverlet << "+-" << stdtverlet << "ms"  << std::endl;


	// write results to file
	for(unsigned int i = 0; i < Neuler.size(); i++) // N.size() == Nana.size()
	filestream << Neuler.at(i) << '\t' // Euler
		<< Nana.at(i) << '\t' // analytical
		<< NeulerMod.at(i) << '\t' // euler modified (runge kutta 2nd order)
		<< Nrk3a.at(i) << '\t' // runge kutta 3rd order
		<< Nverlet.at(i) << '\t' // verlet 
		<< fabs(Neuler.at(i)-Nana.at(i)) << '\t' // error euler
		<< fabs(NeulerMod.at(i)-Nana.at(i)) << '\t' // error mod euler
		<< fabs(Nrk3a.at(i)-Nana.at(i)) << '\t' // error rk
		<< fabs(Nverlet.at(i)-Nana.at(i)) << '\t' // error verlet
		<< std::endl;
	filestream.close();

	//Gnuplot plot
	gp << "set multiplot layout 1,2 \n";
	gp << "plot \"myfile.txt\" using 1 title \" euler \"  with lines, \"myfile.txt\" using 4 title \" RungeKutta3a \" with lines ,  \"myfile.txt\" using 2 title \" analytical\" with lines, \"myfile.txt\" using 5 title \" verlet \" with lines \n ";
	gp << "set logscale y \n";
	gp << "plot \"myfile.txt\" using 7 title \" error euler modified \" with lines, \"myfile.txt\" using 8 title \"error RK3a\" with lines,  \"myfile.txt\" using 6 title \"error Euler\" with lines,  \"myfile.txt\" using 9 title \"error Verlet\" with lines \n";


	return 0;
}

