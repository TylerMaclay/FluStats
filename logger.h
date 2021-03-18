#pragma once

#include <fstream>
#include <string>

//Logger file
//Purpose -- Output stream warnings and events in case of crash or success with strange outcome

class Logger {
private:
	std::ofstream out;

public:
	bool sendMessage(const std::string& message);
	Logger(const std::string& file) : out(file) {}
	~Logger();
};