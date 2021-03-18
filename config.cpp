#include "config.h"

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
					std::string value = dataWithoutComments.substr(equalSign + 1);

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

int determineOutput(std::map<std::string, std::string> configData)
{
	
	int retVal = OutputType::console;
	
	if (configData.find("output") != configData.end()) {
		auto dataVal = configData["output"];
		
		if (dataVal == "csv") {
			retVal = OutputType::csv;
		}
		else if (dataVal == "tsv") {
			retVal = OutputType::tsv;
		}
		else if (dataVal == "txt") {
			retVal = OutputType::txt;
		}
		else if (dataVal == "cout") {
			retVal =  OutputType::console;
		}

	}

	return retVal;
}

std::string determineOutputFile(std::map<std::string, std::string> configData)
{
	if (configData.find("outputFile") != configData.end()) {
		return configData["outputFile"];
	}
	return "./out";
}

std::string determineLogFile(std::map<std::string, std::string> configData)
{
	if (configData.find("logFile") != configData.end()) {
		return configData["logFile"];
	}
	return "./log";
}

double ConfigOptions::getScalingFactor()
{
	return scalingFactor;
}

int ConfigOptions::getStatistic()
{
	return statistic;
}

std::vector<std::string> ConfigOptions::getFileList()
{
	return files;
}

std::string ConfigOptions::getNextFile()
{
	if (fileIt != files.end()) {
		return *fileIt++;
	}
	else return std::string("");
}

int ConfigOptions::getOutput()
{
	return output;
}

std::string ConfigOptions::getOutputFile()
{
	return outputFile;
}

std::string ConfigOptions::getLogfile()
{
	return logFile;
}
