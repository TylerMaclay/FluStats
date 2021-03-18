#include "FS_Stats.h"
#include <cmath>
#include <numeric>

std::vector<double> logAccumulatorFunction(double, int);
std::vector<double> convertLogtoNorm(const std::vector<double>&);

//Generates frequency distribution for MSS-MLE for a given m value
std::vector<double> mssAccumulationFunction(double mGuess, int max = 150) {
	std::vector<double> Accumulator = logAccumulatorFunction(mGuess, max);
	/*Accumulator.push_back(std::expl(-1 * mGuess));
	if (Accumulator[0] == 0) {
		Accumulator[0] = std::numeric_limits<double>::min();
	}
	for (auto i = 1; i < max; i++) {
		double tmp = 0.0;
		for (auto j = 0; j < i; j++) {
			tmp += (Accumulator[j]*mGuess) / static_cast<double>(static_cast<uint64_t>(i) - static_cast<uint64_t>(j) + static_cast<uint64_t>(1));

		}
		tmp /=  i;
		Accumulator.push_back(tmp);
	}*/



	return Accumulator;
}

double logFactorial(double r) {
	if (r == 0 || r == 1) {
		return 0;
	}
	auto retVal = 0.0;
	for (auto i = 1; i < r; i++) {
		retVal += std::log(i);
	}

	return retVal;
}

double addToRollingSum(double sum, double data) {

	const double a = sum;
	const double c = data;
	const double logMaxDouble = std::log(std::numeric_limits<double>::max());
	const double logMinDouble = std::log(std::numeric_limits<double>::min());

	const double logDistance = (c-a);

	if (logDistance > logMaxDouble) {
		return a;
	}
	else if (logDistance < logMinDouble) {
		return c;
	}

	else {
		auto retVal = a + std::log1p(std::exp(logDistance));

		if (std::isinf(retVal) || std::isnan(retVal)) {
			throw std::out_of_range("Sum=" + std::to_string(sum) + "\ndata =" + std::to_string(data) + "\na=" + std::to_string(a) + "\nc=" + std::to_string(c));
		}
		return retVal;
	}
}

std::vector<double> logAccumulatorFunction(double mGuess, int max) {
	std::vector<double> Accumulator;
	
	Accumulator.push_back((-1 * mGuess));
	for (auto i = 1; i < max; i++) {
		double rollingSum = 0.0;
		double constRegion = std::log(mGuess) + Accumulator[0] - std::log(i) - std::log(i+1);
		for (auto j = 1; j < i; j++) {
			double rollingSumTemp = (Accumulator[j] - Accumulator[0] + std::log(i + 1) - std::log(i - j + 1));
			rollingSum = addToRollingSum(rollingSum, rollingSumTemp);
		}

		const double logDoubleMax = std::log(std::numeric_limits<double>::max());
		const double logDoubleMin = std::log(std::numeric_limits<double>::min());

		if (rollingSum > logDoubleMax) {
			Accumulator.push_back(constRegion + rollingSum);
		}
		else if (rollingSum < logDoubleMin) {
			Accumulator.push_back(constRegion);
		}
		else{
			Accumulator.push_back((constRegion + std::log1p(std::exp(rollingSum))));
		}

	}

	return Accumulator;
}

std::vector<double> convertLogtoNorm(const std::vector<double>& old) {
	std::vector<double> newVec;
	for (auto it = old.begin(); it != old.end(); it++) {
		newVec.push_back(std::exp(*it));
	}

	return newVec;
}


//Iteratively looks for maximum likelihood m value. Using a riff on a binary search. Epsilon tells the program which way to go
//Epsilon is required so that program continues to move in the right direction, due to the nature of the distribution generated
//it would be easy to miss the maximal value by averaging current and the greater of high and low.
double mssFindM(double initM, const std::map<int, int>& cultureData) {
	auto maxMutants = (cultureData.rbegin())->first;
	//auto maxVal = (maxMutants > 150 ? (maxMutants + 1) : 150);
	auto maxVal = (maxMutants + 1);
	auto low = 0.00001;
	auto high = static_cast<double>(maxVal);
	auto epsilon = 0.01;
	auto current = maxMutants/std::log(maxMutants);

	auto lowScore = std::async(mssScore,cultureData, mssAccumulationFunction(current - epsilon, maxVal));
	auto currentScore = std::async(mssScore,cultureData, mssAccumulationFunction(current, maxVal));
	auto highScore = std::async(mssScore,cultureData, mssAccumulationFunction(current + epsilon, maxVal));


	lowScore.wait();
	currentScore.wait();
	highScore.wait();


	auto ls = lowScore.get();
	auto cs = currentScore.get();
	auto hs = highScore.get();


	while ((cs < ls || cs < hs)) {
		if (cs < ls) {
			high = current;
			current = (current + low) / 2;
		}
		else if (cs < hs) {
			low = current;
			current = (current + high) / 2;
		}

		lowScore = std::async(mssScore,cultureData, mssAccumulationFunction(current - epsilon, maxVal));
		currentScore = std::async(mssScore,cultureData, mssAccumulationFunction(current, maxVal));
		highScore = std::async(mssScore,cultureData, mssAccumulationFunction(current + epsilon, maxVal));
		
		lowScore.wait();
		currentScore.wait();
		highScore.wait();

		ls = lowScore.get();
		cs = currentScore.get();
		hs = highScore.get();
	}
	return current;
}
std::pair<double, double> getConfidenceInterval(double m, double lstdDev, double cellCount) { //Note! Calculates confidence intervals of the RATE not the mutation count
	double upperCI = std::exp(std::log(m) + (1.96 * lstdDev * std::pow(std::exp(1.96 * lstdDev), -0.315))) / cellCount;
	double lowerCI = std::exp(std::log(m) - (1.96 * lstdDev * std::pow(std::exp(1.96 * lstdDev), +0.315))) / cellCount;
	return std::pair<double, double>(upperCI, lowerCI);
}

double calculateLogStdDev(double m, double cultureCount) {
	return (1.225 * std::pow(m, -0.315) / std::sqrt(cultureCount));
}

std::map<double, double> mSweep(double initM, const std::map<int, int>& cultures) {
	auto low = 750;
	auto high = 780.0;
	auto maxElem = (cultures.rbegin())->first;
	std::map<double, double> ret;
	for (auto i = low; i < high; i++) {
		auto dist = mssAccumulationFunction(i, maxElem +1);
		ret[i] = mssScore(cultures, dist);
	}

	return ret;
}
double mssScore(const std::map<int, int>& equation, const std::vector<double>& distribution) {

	double score = 1.0;
	for (auto it = equation.begin(); it != equation.end(); it++) {
		//auto multScore = distribution[it->first] == 0 ? 1 : distribution[it->first];
		score += (distribution[it->first] * it->second);
	}
	if (std::isinf(score)) {
		score = 0;
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
		auto second = std::stoi(temp.substr(delim + 1));
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

double methodOfTheMedian(double median) {
	return newtonsMethodLCM(median, 1000.0, 0.0, median);
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
	else if (std::abs(upperGuessScore) <= std::abs(lowerGuessScore)) {
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
	return ((r / m) - std::log(m) - 1.24);
}
