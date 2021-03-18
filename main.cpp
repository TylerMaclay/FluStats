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
#include "OutputHandle.h"
#include "logger.h"



int main(int argc, char** argv) {
	try {

		ConfigOptions config(argc, argv, "./config.ini");
		Logger FSLog(config.getLogfile());
		FSLog.sendMessage("Config File Parsed and log file opened!");
		printFileList(config.getFileList());
		std::string outputFile = config.getOutputFile();
		FSLog.sendMessage( "Output file opened " + config.getOutputFile() );
		std::cout << "Scaling Factor: " << config.getScalingFactor() << "\n";
		FSLog.sendMessage("Scaling factor: " + std::to_string(config.getScalingFactor()) );
		std::cout << "Method: " << config.getStatistic() << "\n";
		FSLog.sendMessage("Method: " + std::to_string(config.getStatistic()));
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
			std::string fileProcessMsg = ("Processing: " + file);
			std::cout << fileProcessMsg << std::endl;
			FSLog.sendMessage(fileProcessMsg);
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
			auto cultureCount = 0;
			for (auto it = cultureData.begin(); it != cultureData.end(); it++) {
				cultureCount += it->second;
			}

			auto rate = m / cellCount;
			auto lnstdDev = calculateLogStdDev(m, cultureCount);

			auto confidenceIntervals = getConfidenceInterval(m, lnstdDev, cellCount);
		
			accumulatedData.push_back(OutputContainer(file, m, rate * std::pow(10, config.getScalingFactor()), cellCount, lnstdDev, confidenceIntervals.first * std::pow(10, config.getScalingFactor()), confidenceIntervals.second * std::pow(10, config.getScalingFactor())));
			file = config.getNextFile();
			std::cout << "Done!" << std::endl;
			FSLog.sendMessage("Done!");
		}

		FSLog.sendMessage("Writing data to file...");

		for (auto it = accumulatedData.begin(); it != accumulatedData.end(); it++) {
			output->writeData(*it);
		}

		FSLog.sendMessage("Done!");
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}
