#pragma once
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <future>
#include <limits>

//This file will contain all of the statistics functions for this program to function

class Datapoint {
public:
	int r;
	int N;

	Datapoint(int init_r, int init_N) : r(init_r), N(init_N) {}
};

//LC Method of the median functions and accessories
double methodOfTheMedian(double median);
double newtonsMethodLCM(double initM, double upperBound, double lowerBound, double r);
double leaCoulsonScore(double r, double m);
//End LC Method of the median functions/accessories

//MSS-Maximum Likelihood functions
std::vector<double> mssAccumulationFunction(double mGuess, int max);
std::map<double, double> mSweep(double initM, const std::map<int, int>& cultures);
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

///////Lea-Coulson Method of the Median Model

///////Lea-Coulson Method of the Median Model end

//To Implement: MSS-Maximum Likelihood method
/*
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2041832/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2932672/
https://www.sciencedirect.com/science/article/pii/S0076687905090270?via%3Dihub
*/

