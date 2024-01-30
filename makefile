# some code taken from
# https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# https://stackoverflow.com/questions/313778/generate-dependencies-for-a-makefile-for-a-project-in-c-c

# options for commands
RM = rm -f
CXX = g++
# compile with optimization and no debugging
#CXXFLAGS = -Wall -Werror -Wextra -std=gnu++20 -O3 -static
# compile with debugging and optimization
CXXFLAGS = -Wall -Werror -Wextra -std=gnu++20 -g -O3 -static
# compile with debugging and no optimization
#CXXFLAGS = -Wall -Werror -Wextra -std=gnu++20 -g -static

# flags for boost 1.84 (downloaded by this makefile)
#INCLUDES = -I boost/boost_1_84_0
#LFLAGS = -L boost/boost_1_84_0/stage/lib
#LIBS = -lm -lpng -lboost_program_options -lpthread -lz

# flags for using OS boost version
# ubuntu 22.04 has boost 1.74, need >= 1.80 for C++20
INCLUDES =
LFLAGS =
LIBS = -lm -lpng -lboost_program_options -lpthread -lz

# source files
SRCS_MAIN = $(shell find src/ -maxdepth 1 -name "*.cpp")
SRCS_NONMAIN = $(shell find src/ -mindepth 2 -name "*.cpp")
SRCS_ALL = $(shell find src/ -name "*.cpp")

# dependency and object files
DEPS = $(subst src/,obj/,$(subst .cpp,.d,$(SRCS_ALL)))
DEPS_NONMAIN = $(subst src/,obj/,$(subst .cpp,.d,$(SRCS_NONMAIN)))
OBJS = $(subst src/,obj/,$(subst .cpp,.o,$(SRCS_ALL)))
OBJS_NONMAIN = $(subst src/,obj/,$(subst .cpp,.o,$(SRCS_NONMAIN)))

FFR_BUF = ffr-buf.out
FFR_IMG = ffr-img.out

# buffer renderer executable
$(FFR_BUF): $(OBJS_NONMAIN) obj/ffr_buf.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(FFR_BUF) \
		$(OBJS_NONMAIN) obj/ffr_buf.o $(LFLAGS) $(LIBS)

# image renderer executable
$(FFR_IMG): $(OBJS_NONMAIN) obj/ffr_img.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(FFR_IMG) \
		$(OBJS_NONMAIN) obj/ffr_img.o $(LFLAGS) $(LIBS)

all: $(FFR_BUF) $(FFR_IMG)

# json library
nlohmann/json.hpp:
	mkdir -p nlohmann
	wget -O nlohmann/json.hpp \
	https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp


# need later version of boost because
# the version in Ubuntu 22 does not compile with C++20
#boost:
#	mkdir -p boost/
#	cd boost
#	wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.bz2
#	tar -xvf boost_1_84_0.tar.bz2
#	cd boost_1_84_0/
#	./bootstrap.sh
#	./b2
#	cd ../..
# need to run "sudo ./b2 install --prefix=/usr/local"
# static linking is not a good solution for development
# because glibc shows memory errors in valgrind

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
	$(RM) ./*.out

# so make is not confused if there is a file named "all" or "clean"
.PHONY: all clean
