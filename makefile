# some code taken from
# https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# https://stackoverflow.com/questions/313778/generate-dependencies-for-a-makefile-for-a-project-in-c-c

# options for commands
RM = rm -f
CXX = g++
CXXFLAGS = -Wall -Werror -Wextra -std=gnu++11 -O3 -g
INCLUDES =
LFLAGS =
LIBS = -lm -lpng -lboost_program_options -lpthread

# source files
SRCS_MAIN = $(shell find src/ -maxdepth 1 -name "*.cpp")
SRCS_NONMAIN = $(shell find src/ -mindepth 2 -name "*.cpp")
SRCS_ALL = $(shell find src/ -name "*.cpp")

# dependency and object files
DEPS = $(subst .cpp,.d,$(SRCS_ALL))
OBJS = $(subst .cpp,.o,$(SRCS_ALL))

all: ffbuf.out

ffbuf.out: $(DEPS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o ffbuf.out $(OBJS) $(LFLAGS) $(LIBS)

# json library
nlohmann/json.hpp:
	mkdir -p nlohmann
	wget -O nlohmann/json.hpp \
	https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp

# compile objects
src/%.o: src/%.cpp src/%.d
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

# make dependency files
src/%.d: src/%.cpp nlohmann/json.hpp
	$(CXX) $(CXXFLAGS) -MM -MT $(subst .d,.o,$@) $< -MF $@

include $(DEPS)

clean:
	$(RM) $(OBJS)
	$(RM) $(DEPS)
	$(RM) ffbuf.out

# so make is not confused if there is a file named "all" or "clean"
.PHONY: all clean
