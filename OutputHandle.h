#pragma once

#include <fstream>
#include <string>
#include <iostream>

struct finalData {
	double median;
	double m;
	double i;
};

class OutputContainer {
public:
	std::string fileName;
	double m;
	double rate;
	double cellCount;
	double stdDev;
	double upperCI;
	double lowerCI;

	OutputContainer(std::string name, double _m, double _rate, double _cellCount, double _stdDev, double _upperCI, double _lowerCI) : fileName(name), m(_m), rate(_rate), cellCount(_cellCount), stdDev(_stdDev), upperCI(_upperCI), lowerCI(_lowerCI) {}
};

class OutputHandle {
public:
	virtual void writeData(const OutputContainer& data) =0;
	virtual void writeData(const finalData& data) {};
};

class OutputCSV : public OutputHandle {
	std::ofstream out;
public:
	OutputCSV(const std::string& outFileName);
	void writeData(const OutputContainer& data);
	void writeData(const finalData& data) {
		out << data.i << "," << data.m << "," << data.median << "\n";
	}
};

class OutputTSV : public OutputHandle {
	std::ofstream out;
public:
	OutputTSV(const std::string& outFileName);
	void writeData(const OutputContainer& data);
};

class OutputTXT : public OutputHandle {
	std::ofstream out;
public:
	OutputTXT(const std::string& outFileName);
	void writeData(const OutputContainer& data);
};

class OutputConsole : public OutputHandle {
	void writeData(const OutputContainer& data);
};