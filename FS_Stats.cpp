#include "FS_Stats.h"

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
	auto maxVal = (maxMutants > 150 ? (maxMutants + 1) : 150);
	auto low = 0.00001;
	auto high = static_cast<double>(maxVal);
	auto epsilon = 0.0001;
	auto current = initM;

	auto lowScore = std::async(mssScore,cultureData, mssAccumulationFunction(current - epsilon, maxVal));
	auto currentScore = std::async(mssScore,cultureData, mssAccumulationFunction(current, maxVal));
	auto highScore = std::async(mssScore,cultureData, mssAccumulationFunction(current + epsilon, maxVal));

	lowScore.wait();
	currentScore.wait();
	highScore.wait();

	while (currentScore.get() < lowScore.get() || currentScore.get() < highScore.get()) {
		if (currentScore.get() < lowScore.get()) {
			high = current;
			current = (current + low) / 2;
		}
		else if (currentScore.get() < highScore.get()) {
			low = current;
			current = (current + high) / 2;
		}

		lowScore = std::async(mssScore,cultureData, mssAccumulationFunction(current - epsilon, maxVal));
		currentScore = std::async(mssScore,cultureData, mssAccumulationFunction(current, maxVal));
		highScore = std::async(mssScore,cultureData, mssAccumulationFunction(current + epsilon, maxVal));
		
		lowScore.wait();
		currentScore.wait();
		highScore.wait();
	}
	return current;
}
std::map<double, double> mSweep(double initM, const std::map<int, int>& cultures) {
	auto low = 0.0001;
	auto high = 1000.0;
	auto maxElem = (cultures.rbegin())->first;
	std::map<double, double> ret;
	for (auto i = low; i < high; i += ((high - low) / 100000)) {
		auto dist = mssAccumulationFunction(i, (maxElem > 150 ? (maxElem + 1) : 150));
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
