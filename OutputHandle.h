#pragma once

#include <fstream>
#include <string>
#include <iostream>

class OutputContainer {
public:
	std::string fileName;
	double m;
	double rate;
	double cellCount;

	OutputContainer(std::string name, double _m, double _rate, double _cellCount) : fileName(name), m(_m), rate(_rate), cellCount(_cellCount) {}
};

class OutputHandle {
public:
	virtual void writeData(const OutputContainer& data) =0;
};

class OutputCSV : public OutputHandle {
	std::ofstream out;
public:
	OutputCSV(const std::string& outFileName);
	void writeData(const OutputContainer& data);
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