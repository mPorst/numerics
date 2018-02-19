#include "statistics.h"
#include <vector>
#include <cmath>
#include <chrono>

template <class type>
double mean(std::vector<type> inputArray){
	double temp=0;
	for(unsigned int i=0; i<sizeof(inputArray); i++){
		temp += inputArray.at(i).count();
	}
	return temp/sizeof(inputArray);
}
	
template <class type>
double stdDev(std::vector<type> inputArray, double mean){
	double temp=0;
	for(unsigned int i=0; i<sizeof(inputArray); i++){
		temp += pow(inputArray.at(i).count()-mean,2);
	}
	temp /= (sizeof(inputArray)-1);
	return sqrt(temp);
}


//template double mean<double>(std::vector<double>);
template double mean<std::chrono::duration<double, std::milli>>(std::vector<std::chrono::duration<double, std::milli>>);

//template double stdDev<double>(std::vector<double>, double);
template double stdDev<std::chrono::duration<double, std::milli>>(std::vector<std::chrono::duration<double, std::milli>>, double);
