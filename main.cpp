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
#include <chrono>

//Local Includes
#include "FS_Stats.h"
#include "config.h"
#include "OutputHandle.h"

struct TimeVals {
	std::vector<std::chrono::time_point<std::chrono::steady_clock>> rawVals;
	std::vector<double> adjustedVals;
	std::vector<double> adjacentVals;
	void getTime() {
		if (rawVals.empty()) {
			rawVals.push_back(std::chrono::steady_clock::now());
		}
		else {
			rawVals.push_back(std::chrono::steady_clock::now());
			std::chrono::duration<double> diffAdjust = rawVals.back() - rawVals[0];
			std::chrono::duration<double> diffAdjacent = rawVals.back() - rawVals[rawVals.size() - 2]; //std::vector::size returns number of elements in a vector size-2 gives the second-to-last value
			adjustedVals.push_back(diffAdjust.count());
			adjacentVals.push_back(diffAdjacent.count());
		}
	
	
	};
};

int main(int argc, char** argv) {
	try {
		TimeVals program;
		TimeVals stats;
		program.getTime();
		ConfigOptions config(argc, argv, "./config.ini");
		printFileList(config.getFileList());
		std::string outputFile = config.getOutputFile();
		std::cout << "Scaling Factor: " << config.getScalingFactor() << "\n";
		std::cout << "Method: " << config.getStatistic() << "\n";
		std::vector<OutputContainer> accumulatedData;
		auto file = config.getNextFile();
		OutputHandle* output = nullptr;

		switch (config.getOutput()) {

		case OutputType::csv:
			output = new OutputCSV(outputFile);
			break;
		case OutputType::tsv:
			output = new OutputTSV(outputFile);
			break;
		case OutputType::txt:
			output = new OutputTXT(outputFile);
			break;
		case OutputType::console:
			output = new OutputConsole();
			break;
		default:
			break;
		}

		while (!file.empty()) {
			auto inputData = parseCSV(file);
			auto cultureData = accumulateCultures(inputData);
			auto cellCount = averageCellcounts(inputData);
			auto m = 0.0;
			stats.getTime();
			switch (config.getStatistic()) {
			case Statistics::freq:
				m = findMedian(inputData);
				break;
			case Statistics::lc_mm:
				m = methodOfTheMedian(findMedian(inputData));
				break;
			case Statistics::mss:
				m = mssFindM(methodOfTheMedian(findMedian(inputData)), cultureData);
				break;
			}

			auto rate = m / cellCount;
			accumulatedData.push_back(OutputContainer(file, m, rate * std::pow(10, config.getScalingFactor()), cellCount));
			file = config.getNextFile();
		}

		stats.getTime();

		for (auto it = accumulatedData.begin(); it != accumulatedData.end(); it++) {
			output->writeData(*it);
		}

		for (auto it = stats.adjacentVals.begin(); it != stats.adjacentVals.end(); it++) {
			std::cout << *it << "s\n";
		}
		program.getTime();
		std::cout<<program.adjustedVals[0]<<"s";

		delete output;
		return 0;
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}
