#include "FS_Stats.h"
#include <cmath>
#include <numeric>
#include <iostream>

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
double mssFindM(double initM, const std::map<int, int>& cultureData, int maxIter) {
	auto maxMutants = (cultureData.rbegin())->first;
	auto maxVal = (maxMutants + 1);
	auto currentM = initM;
	auto nextM = initM;
	auto logh = logHsequence(maxVal);
	auto dconvlogh = convolveLogSeries(logh, logh);
	int currIter = 0;
	do {
		currentM = nextM;
		auto currentDist = convertDoubleVecToComplexVec(mssAccumulationFunction(currentM, maxVal));
		std::cout << "Computed distribution, iter: " << currIter << "\n";
		auto firstdiv = std::async(convolveLogSeries<std::complex<double>>, currentDist, logh);
		auto seconddiv = std::async(convolveLogSeries<std::complex<double>>, currentDist, dconvlogh);

		firstdiv.wait();
		seconddiv.wait();

		auto fd = firstdiv.get();
		auto sd = seconddiv.get();

		auto score = scoreFunction(cultureData, currentDist, fd);
		auto info = informationFunction(cultureData, currentDist, fd, sd);
		std::cout << "Score: " << score << "\n";
		std::cout << "Info: " << info << "\n";
		nextM = currentM + (score / info);
		std::cout << "Current M: " << currentM << " nextM: " << nextM << "\n";
		currIter++;
	} while (std::abs(nextM - currentM) > 0.001 && currIter < maxIter);

	return currentM;

}

double scoreFunction(const std::map<int, int>& cultureData, const std::vector<std::complex<double>>& p, const std::vector<std::complex<double>>& dp) {
	auto sum = 0.0;
	
	for (auto it = cultureData.begin(); it != cultureData.end(); it++) {
		sum += std::exp((dp[it->first] - p[it->first])).real();
	}

	return sum;
}

double informationFunction(const std::map<int, int>& cultureData, const std::vector<std::complex<double>>& p, const std::vector<std::complex<double>>& dp, const std::vector<std::complex<double>>& ddp) {
	auto sum = 0.0;
	for (auto it = cultureData.begin(); it != cultureData.end(); it++) {
		sum += std::pow(std::exp(dp[it->first] - p[it->first]).real(), 2) - std::exp(ddp[it->first] - p[it->first]).real();
	}

	return sum;

}

double mssFindMForLikelihoodRatio(const std::map<int, int>& cultureDataOne, const std::map<int, int>& cultureDataTwo, double initM, double cellRatio, int maxIter) {
	auto maxMutantsOne = (cultureDataOne.rbegin())->first;
	auto maxMutantsTwo = (cultureDataTwo.rbegin())->first;
	int currentIter = 0;
	auto maxMutants = maxMutantsOne > maxMutantsTwo ? maxMutantsOne : maxMutantsTwo;
	auto maxVal = (maxMutants + 1);
	auto hseq = logHsequence(maxVal);
	auto dchseq = convolveLogSeries(hseq, hseq);
	auto currentM = initM;
	auto nextM = initM;
	do {
		currentM = nextM;
		auto currentDistOne = convertDoubleVecToComplexVec(mssAccumulationFunction(currentM, maxVal));
		auto currentDistTwo = convertDoubleVecToComplexVec(mssAccumulationFunction(currentM * cellRatio, maxVal));
		std::cout << "Computed distribution, iter: " << currentIter << "\n";
		auto firstDivOne = std::async(convolveLogSeries<std::complex<double>>, currentDistOne, hseq);
		auto firstDivTwo = std::async(convolveLogSeries<std::complex<double>>, currentDistTwo, hseq);
		auto secondDivOne = std::async(convolveLogSeries<std::complex<double>>, currentDistOne, dchseq);
		auto secondDivTwo = std::async(convolveLogSeries<std::complex<double>>, currentDistTwo, dchseq);

		firstDivOne.wait();
		firstDivTwo.wait();
		secondDivOne.wait();
		secondDivTwo.wait();

		auto fdo = firstDivOne.get();
		auto fdt = firstDivTwo.get();
		auto sdo = secondDivOne.get();
		auto sdt = secondDivTwo.get();

		auto score = 0.0;
		auto info = 0.0;

		for (auto it = cultureDataOne.begin(); it != cultureDataOne.end(); it++) {
			score += std::exp((fdo[it->first] - currentDistOne[it->first])).real();
			info += std::pow(std::exp(fdo[it->first] - currentDistOne[it->first]).real(), 2) - std::exp(sdo[it->first] - currentDistOne[it->first]).real();
		}
		for (auto it = cultureDataTwo.begin(); it != cultureDataTwo.end(); it++) {
			score += std::exp((fdt[it->first] - currentDistTwo[it->first])).real();
			info += std::pow(std::exp(fdt[it->first] - currentDistTwo[it->first]).real(), 2) - std::exp(sdt[it->first] - currentDistTwo[it->first]).real();
		}

		std::cout << "Score: " << score << "\n";
		std::cout << "Info: " << info << "\n";

		nextM = currentM + (score / info);
		std::cout << "Current M: " << currentM << " nextM: " << nextM << "\n";
		currentIter++;

	} while (std::abs(currentM - nextM) > 0.001 && currentIter < maxIter);

	currentM = nextM;


	return currentM;
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

	double score = 0.0;
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

double likelihood_ratio(const std::map<int, int>& cultureDataOne, const std::map<int, int>& cultureDataTwo, double mOne, double mTwo, double cellCountOne, double cellCountTwo) {
	
	//Ratio of cell counts to normalize between comparisons
	double cellCountRatio = cellCountOne / cellCountTwo;

	std::cout << "Cell count ratio: " << cellCountRatio << "\n";

	auto likelihood_separate = mssScore(cultureDataOne, mssAccumulationFunction(mOne, ((cultureDataOne.rbegin())->first) + 1)) +
		mssScore(cultureDataTwo, mssAccumulationFunction(mTwo, ((cultureDataTwo.rbegin())->first) + 1));
	std::cout << "Likelihood Separate: " << likelihood_separate << "\n";

	//Calculate null hypothesis 
	//Start by normalizing m2 to m1
	auto combinedM = mssFindMForLikelihoodRatio(cultureDataOne, cultureDataTwo, mOne, cellCountRatio);

	std::cout << "Combined M: " << combinedM << "\n";
	auto nullHypothesis = mssScore(cultureDataOne, mssAccumulationFunction(combinedM, (cultureDataOne.rbegin()->first) + 1)) +
		mssScore(cultureDataTwo, mssAccumulationFunction(combinedM*cellCountRatio, (cultureDataTwo.rbegin()->first) + 1));

	std::cout << "null: " << nullHypothesis << "\n";
	auto likelihood = 2 * (likelihood_separate - nullHypothesis);

	return likelihood; //Outputs the critical value. Need to check against lookup table with 1 dof
	//This is adapted from https://link.springer.com/article/10.1007/s10709-016-9904-3
}

std::vector<std::complex<double>> logHsequence(int max) {
	auto retval = std::vector<std::complex<double>>{};
	retval.push_back(std::log(std::complex<double>(-1)));
	for (auto it = 1; it < max; it++) {
		auto diff = (1.0 / it) - (1.0 / (it + 1));
		retval.push_back(std::log(std::complex<double>(diff, 0)));
	}

	return retval;
}

std::vector<std::complex<double>> convertDoubleVecToComplexVec(const std::vector<double>& v) {
	auto retval = std::vector<std::complex<double>>{};
	for (auto i : v) {
		retval.push_back(std::complex<double>(i));
	}

	return retval;
}


double incompleteGamma(double s, double x) {
	auto retval = (s * std::log(x)) - x + std::lgamma(s);
	int k = 0;
	auto testval = (k * std::log(x) - std::lgamma(s + k + 1));
	auto rollingSum = std::exp(testval);
	while (testval > -10) {
		k++;
		testval = k * std::log(x) - std::lgamma(s + k + 1);
		rollingSum += std::exp(testval);
	}

	retval += std::log(rollingSum);

	return std::exp(retval);
}


double chiSquareCDF(double x, double dof) {
	return (incompleteGamma(dof / 2.0, x / 2.0) / tgamma(dof / 2.0));
}