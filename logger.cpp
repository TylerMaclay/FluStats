#include "logger.h"

Logger::~Logger() {
	if(out.is_open()){
		out.close();
	}
}

bool Logger::sendMessage(const std::string& message) {
	out << message << std::endl;
	out.flush();
	return out.fail();
}