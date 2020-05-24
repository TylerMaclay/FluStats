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



//Config file input
std::map<std::string, std::string> parseConfigFile(const std::string& file);
void printConfig(std::map<std::string, std::string> data);
std::vector<std::string> fileListReader(int argc, char** argv);
void printFileList(std::vector<std::string> fileList);
double determineScalingFactor(std::map<std::string, std::string> configData);
int determineStatistic(std::map<std::string, std::string> configData);

int main(int argc, char** argv) {
	try {
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
		}

		return 0;
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}



//Science and statistics stop here, below is program administration
std::map<std::string, std::string> parseConfigFile(const std::string& file) {
	std::ifstream in(file);
	std::map<std::string, std::string> parsedValues;
	if (!in.is_open()) {
		throw std::runtime_error("Could not open config file, check that file is valid!");
	}
	else {
		std::string rawData = "";
		while (std::getline(in, rawData)) {
			//Comments begin with ";", this subsets away from that.
			std::string dataWithoutComments = rawData.substr(0, rawData.find(';'));
			//Removes all whitespace from string then checks if the string is empty. Skipping the line if it is:
			dataWithoutComments.erase(std::remove_if(dataWithoutComments.begin(), dataWithoutComments.end(), isspace), dataWithoutComments.end());
			if (dataWithoutComments.empty()) {
				continue;
			}
			else {
				auto equalSign = dataWithoutComments.find('=');
				if (equalSign == std::string::npos) {
					throw std::runtime_error("Improperly formatted config file. Lines should be key=value");
				}
				else {
					std::string key = dataWithoutComments.substr(0, equalSign);
					std::string value = dataWithoutComments.substr(equalSign+1);

					if (parsedValues.find(key) != parsedValues.end()) {
						std::cerr << "Warning: a setting has been defined twice. Using most recent definition";
					}

					parsedValues[key] = value;
				}

			}

		}
	}

	return parsedValues;
}
void printConfig(std::map<std::string, std::string> data) {
	for (auto it = data.begin(); it != data.end(); it++) {
		std::cout << it->first << "\t" << it->second << "\n";
	}
}

std::vector<std::string> fileListReader(int argc, char** argv)
{
	std::vector<std::string> data;
	for (auto i = 1; i < argc; i++) {
		data.push_back(std::string(argv[i]));
	}

	return data;
}

void printFileList(std::vector<std::string> fileList)
{
	for (auto i : fileList) {
		std::cout << i << "\n";
	}
}

double determineScalingFactor(std::map<std::string, std::string> configData)
{
	//If a scaling factor is defined, use that if not default to x10^6
	if (configData.find("scale") != configData.end()) {
		return std::stod((configData.find("scale"))->second);
	}
	return 6.0;
}

int determineStatistic(std::map<std::string, std::string> configData) //Using enum here to make code more readable. Translates to mss == 0, lc_mm == 1, freq == 2
{
	if (configData.find("method") != configData.end()) {
		auto methodVal = (configData.find("method"))->second;
		if (methodVal == "mss") {
			return Statistics::mss;
		}
		else if (methodVal == "lc-mm") {
			return Statistics::lc_mm;
		}
		else if (methodVal == "freq") {
			return Statistics::freq;
		}
		else {
			throw std::runtime_error("Invalid statistic selection! Please choose one of (mss) (lc-mm) (freq)");
		}
	}
	else {
		return Statistics::mss;
	}
}
