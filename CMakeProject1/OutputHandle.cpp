#include "OutputHandle.h"

OutputCSV::OutputCSV(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error! Could not open output file!");
	}
	out << "Input Filename," << "Mutation value," << "Cell Count," << "Rate," << "StdDev - ln(m)," << "Upper 95% CL," << "Lower 95% CL" <<"\n";
}

void OutputCSV::writeData(const OutputContainer& data)
{
	out << data.fileName << "," << data.m << "," << data.cellCount << "," << data.rate << "," << data.stdDev << "," << data.upperCI << "," << data.lowerCI << "\n";
}

OutputTSV::OutputTSV(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error! Cannot Open File");
	}
	out << "Input Filename\t" << "Mutation value\t" << "Cell Count\t" << "Rate\t" << "StdDev - ln(m)\t" << "Upper 95% CL\t" << "Lower 95% CL\n";
}

void OutputTSV::writeData(const OutputContainer& data)
{
	out << data.fileName << "\t" << data.m << "\t" << data.cellCount << "\t" << data.rate << "\t" << data.stdDev << "\t" << data.upperCI << "\t" << data.lowerCI << "\n";
}

OutputTXT::OutputTXT(const std::string& outFileName)
{
	out.open(outFileName);
	if (!out.is_open()) {
		throw std::runtime_error("Error, could not open output file");
	}
	out << "Input Filename " << "Mutation value " << "Cell Count " << "Rate " << "StdDev - ln(m) " << "Upper 95% CL " << "Lower 95% CL\n";
}

void OutputTXT::writeData(const OutputContainer& data)
{
	out << data.fileName << " " << data.m << " " << data.cellCount << " " << data.rate << " " << data.stdDev << " " << data.upperCI << " " << data.lowerCI << "\n";
}

void OutputConsole::writeData(const OutputContainer& data)
{
	std::cout << data.fileName << ": " << data.m << " " << data.cellCount << " " << data.rate << " " << data.stdDev << " " << data.upperCI << " " << data.lowerCI << "\n";
}
