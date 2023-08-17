#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>


//Makes statistics more readable than ints
enum Statistics { mss, lc_mm, freq };
enum OutputType { csv, tsv, txt, console };

//Config file input
std::map<std::string, std::string> parseConfigFile(const std::string& file);
void printConfig(std::map<std::string, std::string> data);
std::vector<std::string> fileListReader(int argc, char** argv);
void printFileList(std::vector<std::string> fileList);
double determineScalingFactor(std::map<std::string, std::string> configData);
int determineStatistic(std::map<std::string, std::string> configData);
int determineOutput(std::map<std::string, std::string> configData);
std::string determineOutputFile(std::map<std::string, std::string> configData);
bool determineCalculateStats(std::map<std::string, std::string> configData);

//Class to hold config options
class ConfigOptions {
private:
	double scalingFactor;
	int statistic;
	int output;
	std::vector<std::string> files;
	std::map<std::string, std::string> configDataRaw;
	std::string configFileLocation;
	std::vector<std::string>::iterator fileIt;
	std::string outputFile;
	bool CalcStats;



public:

	ConfigOptions(int argc, char** argv, const std::string& configFile ="") : configFileLocation(configFile), scalingFactor(1000000.0), statistic(0), output(OutputType::console), outputFile("./out") {
		if (!configFile.empty()) {
			configDataRaw = parseConfigFile(configFile);
			scalingFactor = determineScalingFactor(configDataRaw);
			statistic = determineStatistic(configDataRaw);
			output = determineOutput(configDataRaw);
			outputFile = determineOutputFile(configDataRaw);
		}
		files = fileListReader(argc, argv);
		fileIt = files.begin();

	}

	double getScalingFactor();
	int getStatistic();
	std::vector<std::string> getFileList();
	std::string getNextFile();
	int getOutput();
	std::string getOutputFile();
	bool getCalcStats();
	int getNumFiles();


};


