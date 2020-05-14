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

class Datapoint {
public:
	int r;
	int N;

	Datapoint(int init_r, int init_N) : r(init_r), N(init_N) {}
};
enum Statistics { mss, lc_mm, freq};

//LC Method of the median functions and accessories
double methodOfTheMedian(double median);
double newtonsMethodLCM(double initM, double upperBound, double lowerBound, double r);
double leaCoulsonScore(double r, double m);
//End LC Method of the median functions/accessories

//MSS-Maximum Likelihood functions
std::vector<double> mssAccumulationFunction(double mGuess, int max);
std::map<double, double> mSweep(double initM, const std::map<int,int>& cultures);
double mssScore(const std::map<int, int>& equation, const std::vector<double>& distribution);
double mssFindM(double initM, const std::map<int, int>& cultureData);
//MSS-Maximum Likelihood functions end


///Data input processing
std::vector<Datapoint> parseCSV(const std::string& fileName);
std::map<int, int> accumulateCultures(const std::vector<Datapoint>& data);
//End Data input processing

//Basic stats
double averageCellcounts(const std::vector<Datapoint>& data);
double findMedian(const std::vector<Datapoint>& input);
//End Basic stats

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
		/*
		auto input = parseCSV("./rad59.csv");
		auto cultureData = accumulateCultures(input);
		auto cells = averageCellcounts(input);
		auto m = mssFindM(methodOfTheMedian(findMedian(input)), cultureData);
		auto rate = (m / cells);
		auto tenminussixrate = rate * std::pow(10, 6);
		auto sigma = 1.225 * std::pow(m, -0.315) / std::sqrt(input.size());
		std::cout << "m value: " << mssFindM(methodOfTheMedian(findMedian(input)), cultureData) << "\n";
		std::cout << "Rate: " << tenminussixrate << "x10^-6\n";*/

		return 0;
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}

///////Lea-Coulson Method of the Median Model
double methodOfTheMedian(double median) {
	return newtonsMethodLCM(7, 1000.0, 0.0, median);
}

double newtonsMethodLCM(double initM, double upperBound, double lowerBound, double r) {

	auto upperGuessM = (initM + upperBound) / 2;
	auto lowerGuessM = (initM + lowerBound) / 2;

	auto upperGuessScore = leaCoulsonScore(r, upperGuessM);
	auto lowerGuessScore = leaCoulsonScore(r, lowerGuessM);
	auto currentGuessScore = leaCoulsonScore(r, initM);

	/* Debugging code

	std::cout << "Current Guess: " << initM << "\tCurrent Guess Score: " << currentGuessScore << "\n";
	std::cout << "Lower Guess: " << lowerGuessM << "\tLower Guess Score: " << lowerGuessScore << "\n";
	std::cout << "Upper Guess: " << upperGuessM << "\tUpper Guess Score: " << upperGuessScore << "\n";
	
	*/

	if (initM == upperBound || initM == lowerBound) {
		return initM;
	}
	
	if (std::abs(currentGuessScore) < 0.001) {
		return initM;
	}
	else if (std::abs(upperGuessScore) <= std::abs(lowerGuessScore)){
		return newtonsMethodLCM(upperGuessM, upperBound, lowerGuessM, r);
	}
	else {
		return newtonsMethodLCM(lowerGuessM, upperGuessM, lowerBound, r);
	}
}

//leaCoulsonMedCalc calculates the relationship between number of mutants (r) and mutation rate (m).
double leaCoulsonScore(double r, double m) {
	if (m <= 0) {
		throw std::out_of_range("Potential mutation rate must be >1");
	}
	return ((r/m)-std::log(m)-1.24);
}

///////Lea-Coulson Method of the Median Model end

//To Implement: MSS-Maximum Likelihood method
/*
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2041832/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2932672/
https://www.sciencedirect.com/science/article/pii/S0076687905090270?via%3Dihub
*/

//Generates frequency distribution for MSS-MLE for a given m value
std::vector<double> mssAccumulationFunction(double mGuess, int max = 150) {
	std::vector<double> Accumulator;
	Accumulator.push_back(std::exp(-1 * mGuess));
	for (auto i = 1; i < max; i++) {
		double tmp = 0.0;
		for (auto j = 0; j < i; j++) {
			tmp += Accumulator[j] / static_cast<double>(static_cast<uint64_t>(i) - static_cast<uint64_t>(j) + static_cast<uint64_t>(1));
		}
		tmp *= mGuess / i;
		Accumulator.push_back(tmp);
	}

	return Accumulator;
}

//Iteratively looks for maximum likelihood m value. Using a riff on a binary search. Epsilon tells the program which way to go
//Epsilon is required so that program continues to move in the right direction, due to the nature of the distribution generated
//it would be easy to miss the maximal value by averaging current and the greater of high and low.
double mssFindM(double initM, const std::map<int, int>& cultureData) {
	auto maxMutants = (cultureData.rbegin())->first;
	auto maxVal = (maxMutants > 150 ? (maxMutants+1) : 150);
	auto low = 0.00001;
	auto high = 100.0;
	auto epsilon = 0.001;
	auto current = initM;
	auto lowScore = mssScore(cultureData, mssAccumulationFunction(current - epsilon, maxVal));
	auto currentScore = mssScore(cultureData, mssAccumulationFunction(current, maxVal));
	auto highScore = mssScore(cultureData, mssAccumulationFunction(current + epsilon, maxVal));

	while (currentScore < lowScore || currentScore < highScore) {
		if (currentScore < lowScore) {
			high = current;
			current = (current + low) / 2;
		}
		else if (currentScore < highScore) {
			low = current;
			current = (current + high) / 2;
		}

		lowScore = mssScore(cultureData, mssAccumulationFunction(current - epsilon, maxVal));
		currentScore = mssScore(cultureData, mssAccumulationFunction(current, maxVal));
		highScore = mssScore(cultureData, mssAccumulationFunction(current + epsilon, maxVal));
	}
	return current;
}
std::map<double, double> mSweep(double initM, const std::map<int, int>& cultures) { 
	auto low = 0.0001;
	auto high = initM * 10.0;
	auto maxElem = (cultures.rbegin())->first;
	std::map<double, double> ret;
	for (auto i = low; i < high; i += ((high - low) / 10000)) {
		auto dist = mssAccumulationFunction(i, (maxElem > 150 ? maxElem : 150));
		ret[i] = mssScore(cultures, dist);
	}

	return ret;
}
double mssScore(const std::map<int, int>& equation, const std::vector<double>& distribution) {
	
	double score = 1.0;
	for (auto it = equation.begin(); it != equation.end(); it++) {
		score *= std::pow(distribution[it->first], it->second);
	}

	return score;
}

//Collects data from well formatted csv files. Expected files to have general format of r,N\n
std::vector<Datapoint> parseCSV(const std::string& fileName) {
	std::ifstream in(fileName);
	if (!in.is_open()) {
		throw std::runtime_error("Error! File could not be opened, check filename");
	}
	std::string temp;
	std::vector<Datapoint> data;
	while (std::getline(in, temp)) {
		//Find tab
		auto delim = temp.find(',');
		if (delim == std::string::npos) {
			throw std::logic_error("Error, data file improperly formatted");
		}
		auto first = std::stoi(temp.substr(0, delim));
		auto second = std::stoi(temp.substr(delim+1));
		data.push_back(Datapoint(first, second));
	}

	return data;
}

//Generates histogram data of culture data (maybe should rename function?)
std::map<int, int> accumulateCultures(const std::vector<Datapoint>& data) {
	std::map<int, int> dataAccumulate;
	for (auto it = data.begin(); it != data.end(); it++) {
		auto mutVal = it->r;
		if (dataAccumulate.find(mutVal) == dataAccumulate.end()) {
			dataAccumulate[mutVal] = 1;
		}
		else {
			dataAccumulate[mutVal] += 1;
		}
	}

	return dataAccumulate;
}

//Takes cell counts reported (N) and averages them.
double averageCellcounts(const std::vector<Datapoint>& data) {
	double rollingAverage = 0.0;
	for (auto it = data.begin(); it != data.end(); it++) {
		rollingAverage += (static_cast<double>(it->N) / data.size());
	}

	return rollingAverage;
}
//Finds the median of raw r values. Note! Datapoint is the input here not the histogram generated above.
double findMedian(const std::vector<Datapoint>& input) {
	std::vector<int> cultures;
	for (auto it = input.begin(); it != input.end(); it++) {
		cultures.push_back(it->r);
	}

	std::sort(cultures.begin(), cultures.end());

	auto cultureSize = cultures.size();

	double median = 0;

	if (cultureSize % 2 == 0) {
		median = ((static_cast<double>(cultures[cultureSize / 2]) + static_cast<double>(cultures[(cultureSize / 2) - 1])) / 2.0);
	}
	else {
		median = cultures[(cultureSize + 1) / 2];
	}

	return median;
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
