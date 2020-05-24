#include "OutputHandle.h"

OutputCSV::OutputCSV(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error! Could not open output file!");
	}
	out << "Input Filename," << "Mutation value," << "Cell Count," << "Rate\n";
}

void OutputCSV::writeData(const OutputContainer& data)
{
	out << data.fileName << "," << data.m << "," << data.cellCount << "," << data.rate << "\n";
}

OutputTSV::OutputTSV(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error! Cannot Open File");
	}
	out << "Input Filename\t" << "Mutation value\t" << "Cell Count\t" << "Rate\n";
}

void OutputTSV::writeData(const OutputContainer& data)
{
	out << data.fileName << "\t" << data.m << "\t" << data.cellCount << "\t" << data.rate << "\n";
}

OutputTXT::OutputTXT(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error, could not open output file");
	}
	out << "Input Filename " << "Mutation value " << "Cell Count " << "Rate\n";
}

void OutputTXT::writeData(const OutputContainer& data)
{
	out << data.fileName << " " << data.m << " " << data.cellCount << " " << data.rate << "\n";
}

void OutputConsole::writeData(const OutputContainer& data)
{
	std::cout << data.fileName << ": " << data.m << " " << data.cellCount << " " << data.rate << "\n";
}
