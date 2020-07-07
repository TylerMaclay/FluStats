//Fluctuation analysis tool
//Tyler Maclay
//Algorithms implemented as described in Rosche & Foster. Methods 2000
//Trying to clone: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687991/



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
#include "OutputHandle.h"



int main(int argc, char** argv) {
	try {

		ConfigOptions config(argc, argv, "./config.ini");
		printFileList(config.getFileList());
		std::string outputFile = config.getOutputFile();
		std::cout << "Scaling Factor: " << config.getScalingFactor() << "\n";
		std::cout << "Method: " << config.getStatistic() << "\n";
		std::vector<OutputContainer> accumulatedData;
		auto file = config.getNextFile();
		OutputHandle* output = nullptr;

		//Creates output object
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

		//Below processes files and sends data for appropriate statistics
		while (!file.empty()) {
			auto inputData = parseCSV(file);
			auto cultureData = accumulateCultures(inputData);
			auto cellCount = averageCellcounts(inputData);
			auto m = 0.0;

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

		//Below writes output to whatever was chosen in config
		for (auto it = accumulatedData.begin(); it != accumulatedData.end(); it++) {
			output->writeData(*it);
		}

		delete output;
		return 0;
	}

	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}
