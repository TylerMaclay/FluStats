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



int main(int argc, char** argv) {
	try {

		/*ConfigOptions config(argc, argv, "./config.ini");
		//printFileList(config.getFileList());
		std::string outputFile = config.getOutputFile();
		std::cout << "Scaling Factor: " << config.getScalingFactor() << "\n";
		std::cout << "Method: " << config.getStatistic() << "\n";
		std::vector<OutputContainer> accumulatedData;
		//auto file = config.getNextFile();
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
		

		std::vector<finalData> findat;

		for (auto i = 1.0; i < 20000; i*=1.1) {
			auto cellCount_init = 1E6;

			if (std::find_if(findat.begin(), findat.end(), [&](const finalData& data){return std::round(i) == std::round(data.i); }) != findat.end()) {
				continue;
			}
			std::vector<Datapoint> inputData;
			for (int j = 0; j < 10; j++) {
				inputData.push_back(Datapoint(std::round(i), cellCount_init));
			}
			auto cultureData = accumulateCultures(inputData);
			auto cellCount = averageCellcounts(inputData);
			auto m = 0.0;


			auto median = methodOfTheMedian(findMedian(inputData));
			m = mssFindM(methodOfTheMedian(findMedian(inputData)), cultureData);
			auto cultureCount = 0;
			for (auto it = cultureData.begin(); it != cultureData.end(); it++) {
				cultureCount += it->second;
			}

			findat.push_back(finalData{ median, m, std::round(i) });
		}

		for (auto it = findat.begin(); it != findat.end(); it++) {
			output->writeData(*it);
		}
		*/
		std::vector<Datapoint> ctrl_th;
		std::vector<Datapoint> thirteennn;
		std::vector<Datapoint> fiftyeightff;
		std::vector<Datapoint> oneD;

		for (auto i = 0; i < 10; i++) {
			ctrl_th.push_back(Datapoint(252, 1E6));
			thirteennn.push_back(Datapoint(1399, 1E6));
			fiftyeightff.push_back(Datapoint(5844, 1E6));
			oneD.push_back(Datapoint(1, 1E6));
		}

		auto one = accumulateCultures(ctrl_th);
		auto two = accumulateCultures(thirteennn);
		auto three = accumulateCultures(fiftyeightff);
		auto oneD_t = accumulateCultures(oneD);

		/*auto msone = mSweep(10, one);
		
		std::ofstream outone("out1.csv");
		outone << "m,probability\n";
		for (auto it = msone.begin(); it != msone.end(); it++) {
			outone << it->first << "," << it->second << "\n";
		}

		outone.close();
		*/
		auto mstwo = mSweep(10, two);

		std::ofstream outtwo("out2.csv");
		outtwo << "m,probability\n";
		for (auto it = mstwo.begin(); it != mstwo.end(); it++) {
			outtwo << it->first << "," << it->second << "\n";
		}

		outtwo.close();


		auto msthree = mSweep(10, three);

		std::ofstream outth("out3.csv");
		outtwo << "m,probability\n";
		for (auto it = msthree.begin(); it != msthree.end(); it++) {
			outth << it->first << "," << it->second << "\n";
		}

		outth.close();
		/*
		auto msOned = mSweep(10, oneD_t);
		std::ofstream outd("out4.csv");
		outone << "m,probability\n";
		for (auto it = msOned.begin(); it != msOned.end(); it++) {
			outth << it->first << "," << it->second << "\n";
		}

		outd.close();
		//delete output;

		*/
		return 0;
	}
	catch (const std::exception& e) {
		std::cerr<<e.what();
		return -1;
	}
}
