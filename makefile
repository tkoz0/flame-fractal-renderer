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
DEPS = $(subst src/,obj/,$(subst .cpp,.d,$(SRCS_ALL)))
OBJS = $(subst src/,obj/,$(subst .cpp,.o,$(SRCS_ALL)))
OBJS_NONMAIN = $(subst src/,obj/,$(subst .cpp,.o,$(SRCS_NONMAIN)))

# executable names
FFR_BASIC = ffr-basic.out
FFR_BUFFER = ffr-buffer.out

all: $(FFR_BASIC) $(FFR_BUFFER)

$(FFR_BASIC): $(DEPS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(FFR_BASIC) \
		$(OBJS_NONMAIN) obj/ffr_basic.o $(LFLAGS) $(LIBS)

$(FFR_BUFFER): $(DEPS) $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(FFR_BUFFER) \
		$(OBJS_NONMAIN) obj/ffr_buffer.o $(LFLAGS) $(LIBS)

# json library
nlohmann/json.hpp:
	mkdir -p nlohmann
	wget -O nlohmann/json.hpp \
	https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp

# compile objects
obj/%.o: src/%.cpp obj/%.d
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

# make dependency files
obj/%.d: src/%.cpp .prereq
	$(CXX) $(CXXFLAGS) -MM -MT $(subst src/,obj/,$(subst .d,.o,$@)) $< -MF $@

include $(DEPS)

# prerequisites for entire project
.prereq: nlohmann/json.hpp
	mkdir -p $(shell dirname $(OBJS))
	touch .prereq

clean:
	$(RM) .prereq
	$(RM) $(OBJS)
	$(RM) $(DEPS)
	$(RM) *.out

# so make is not confused if there is a file named "all" or "clean"
.PHONY: all clean
