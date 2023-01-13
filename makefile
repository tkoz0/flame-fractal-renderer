# basic akefile to compile the entire program at once

LIBS=-lm -lpng -lboost_program_options -lpthread

all: nlohmann/json.hpp
	mkdir -p bin
	g++ -Wall -Werror -Wextra -std=gnu++11 -O3 -o bin/ffbuf src/*.cpp $(LIBS)

nlohmann/json.hpp:
	mkdir -p nlohmann
	wget -O nlohmann/json.hpp \
	https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp
