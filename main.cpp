//Fluctuation analysis tool
//Tyler Maclay
//Algorithms implemented as described in Rosche & Foster. Methods 2000
//Trying to clone: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687991/

//ToDo: Create function to calculate chosen statistic across submitted files
//ToDo: Generate an output function so for when more than one file is submitted
//ToDo: Refactor this mess into something readable, maybe alias some commonly used and verbose types


#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>

//Local Includes
#include "FS_Stats.h"
#include "config.h"

int main(int argc, char** argv) {
	try {
		
		ConfigOptions config(argc, argv, "./config.ini");
		printFileList(config.getFileList());
		std::cout << "Scaling Factor: " << config.getScalingFactor() << "\n";
		std::cout << "Method: " << config.getStatistic() << "\n";
		/*
		auto test = parseConfigFile("./config.ini");
		printConfig(test);
		auto fileList = fileListReader(argc, argv);
		printFileList(fileList);
		std::cout << "Scaling Factor: " << determineScalingFactor(test) << "\n";
		std::cout << "Method: " << determineStatistic(test) << "\n";
		
		auto input = parseCSV("./rad59.csv");
		auto cultureData = accumulateCultures(input);
		auto cells = averageCellcounts(input);
		auto m = mssFindM(methodOfTheMedian(findMedian(input)), cultureData);
		auto rate = (m / cells);
		auto tenminussixrate = rate * std::pow(10, 6);
		auto sigma = 1.225 * std::pow(m, -0.315) / std::sqrt(input.size());
		std::cout << "m value: " << mssFindM(methodOfTheMedian(findMedian(input)), cultureData) << "\n";
		std::cout << "Rate: " << tenminussixrate << "x10^-6\n";

		auto rawData = mSweep(findMedian(input), cultureData);
		std::ofstream out("dataOut.csv");
		if (!out.is_open()) { throw std::runtime_error("Could not open output file!"); }
		for (auto it = rawData.begin(); it != rawData.end(); it++) {
			out << it->first << "," << it->second << "\n";
		}*/

		return 0;
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}
