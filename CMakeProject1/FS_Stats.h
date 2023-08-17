#pragma once
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <future>
#include <limits>
#include <complex>

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
double mssFindM(double initM, const std::map<int, int>& cultureData, int maxIter = 20);
std::pair<double, double> getConfidenceInterval(double m, double lstdDev, double cellCount);
double calculateLogStdDev(double m, double cultureCount);

double likelihood_ratio(const std::map<int,int>& cultureDataOne, const std::map<int,int>& cultureDataTwo, double mOne, double mTwo, double cellCountOne, double cellCountTwo);
double mssFindMForLikelihoodRatio(const std::map<int, int>& cultureDataOne, const std::map<int, int>& cultureDataTwo, double initM, double cellRatio, int maxIter= 20);

std::vector<std::complex<double>> logHsequence(int max);
std::vector<std::complex<double>> convertDoubleVecToComplexVec(const std::vector<double>& v);

double scoreFunction(const std::map<int, int>& cultureData, const std::vector<std::complex<double>>& p, const std::vector<std::complex<double>>& dp);
double informationFunction(const std::map<int, int>& cultureData, const std::vector<std::complex<double>>& p, const std::vector<std::complex<double>>& dp, const std::vector<std::complex<double>>& ddp);
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


//Chisq Implementation
double chiSquareCDF(double x, double dof);
double incompleteGamma(double s, double x);

template <typename T>
std::vector<T> convolveLogSeries(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    auto returnVal = std::vector<T>{};
    if (rhs.size() == lhs.size()) {
        for (auto n = 0; n < rhs.size(); n++) {
            auto sum = rhs[0] + lhs[n];
            for (auto i = 1; i <= n; i++) {
                auto tempsum = rhs[i] + lhs[n - i] - sum;
                if (std::abs(tempsum) > std::log(std::numeric_limits<double>::min()) && std::abs(tempsum) < std::log(std::numeric_limits<double>::max())) {
                    sum += std::log(1.0 + std::exp(tempsum)); //does not add to sum if number is smaller than lower limit of resolution
                }
                else if (std::abs(tempsum) < std::log(std::numeric_limits<double>::min())) {
                    continue;
                }
                else if (std::abs(tempsum) > std::log(std::numeric_limits<double>::max())) {
                    sum += tempsum;
                }
            }

            returnVal.push_back(sum);
        }
    }

    return returnVal;
}