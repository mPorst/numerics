#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "gnuplot-iostream.h"
#include "numerics.h"
#include <chrono>
#include "statistics.h"

int main(){
	
	//set up gnuplot interface
	Gnuplot gp;

	double stepsize = 500;
	double numberMolecules = 100000000;
	double tau1 = 100000;
	std::vector<double> Nana, Neuler, NeulerMod, Nrk3a;
	
	std::vector<std::chrono::duration<double, std::milli>> teuler, tana, teulerMod, trk3a;
	std::fstream filestream;
	filestream.open("myfile.txt", std::fstream::out);

	// calculate 1000 times to account for variance of time measurements
	for(int i=0; i<1000; i++) {
	// call numerical methods and measure time points between them
	//
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	Neuler = euler(stepsize, numberMolecules, tau1, 1);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	Nana = analytical(stepsize, numberMolecules, tau1, 1);
	std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
	NeulerMod = eulerModified(stepsize, numberMolecules, tau1, 1);
	std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
	Nrk3a = RK3a(stepsize, numberMolecules, tau1, 1);
	std::chrono::high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
	
	//time point to duration cast
	teuler.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1));
	tana.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2));
	teulerMod.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3));
	trk3a.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4));
	}

	//mean value and standard deviation
	
	// get mean of execution time and standard deviation (defined in statistics)
	double meanteuler=0, meantana=0, meanteulerMod=0, meantrk3a=0;
	double stdteuler=0, stdtana=0, stdteulerMod=0, stdtrk3a=0;
	
	meanteuler = mean<std::chrono::duration<double,std::milli>>(teuler);
	stdteuler = stdDev<std::chrono::duration<double,std::milli>>(teuler, meanteuler);
	meantana = mean<std::chrono::duration<double,std::milli>>(tana);
	stdtana = stdDev<std::chrono::duration<double,std::milli>>(tana, meantana);
	meanteulerMod = mean<std::chrono::duration<double,std::milli>>(teulerMod);
	stdteulerMod = stdDev<std::chrono::duration<double,std::milli>>(teulerMod, meanteulerMod);
	meantrk3a = mean<std::chrono::duration<double,std::milli>>(trk3a);
	stdtrk3a = stdDev<std::chrono::duration<double,std::milli>>(trk3a, meantrk3a);

	std::cout << "Euler: " << meanteuler << "+-" << stdteuler << "ms"  << std::endl
		<< "analytical: " << meantana << "+-" << stdtana  << "ms" << std::endl
		<< "Runge Kutta 2nd order(mod. Euler): " << meanteulerMod << "+-" << stdteulerMod << "ms"  << std::endl
		<< "Runge Kutta 3rd order: " << meantrk3a << "+-" << stdtrk3a << "ms"  << std::endl;


	// write results to file
	for(unsigned int i = 0; i < Neuler.size(); i++) // N.size() == Nana.size()
	filestream << Neuler.at(i) << '\t' // Euler
		<< Nana.at(i) << '\t' // analytical
		<< NeulerMod.at(i) << '\t' // euler modified (runge kutta 2nd order)
		<< Nrk3a.at(i) << '\t' // runge kutta 3rd order
		<< fabs(Neuler.at(i)-Nana.at(i)) << '\t' // error euler
		<< fabs(NeulerMod.at(i)-Nana.at(i)) << '\t' // error mod euler
		<< fabs(Nrk3a.at(i)-Nana.at(i)) << '\t' // error rk
		<< std::endl;
	filestream.close();

	//Gnuplot plot
	gp << "set multiplot layout 1,2 \n";
	gp << "plot \"myfile.txt\" using 3 title \" euler \"  with lines, \"myfile.txt\" using 4 title \" RungeKutta3a \" with lines ,  \"myfile.txt\" using 2 title \" analytical\" with lines \n ";
	gp << "set logscale y \n";
	gp << "plot \"myfile.txt\" using 6 title \" error euler modified \" with lines, \"myfile.txt\" using 7 title \"error RK3a\" with lines,  \"myfile.txt\" using 5 title \"error Euler\" with lines \n";


	return 0;
}

