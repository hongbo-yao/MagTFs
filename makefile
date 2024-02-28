CXX					:=  g++
CXXFLAGS				:=  -std=c++11 -O3 -fopenmp

INCLUDES				:= $(wildcard src/*.h)
SRCS					:= $(wildcard src/*.cpp) 
OBJS					:= $(patsubst %.cpp, %.o, $(SRCS))

MagTFs_INCLUDE		:= -I src/
EIGEN_INCLUDE		:= -I contrib/eigen

%.o : %.cpp 
	@echo "MagTFs is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(MagTFs_INCLUDE) $(EIGEN_INCLUDE) -c $< -o $@

MagTFs: $(OBJS) $(INCLUDES)
	@$(CXX) $(CXXFLAGS) -o MagTFs $(OBJS)
	
clean:	
	@rm -rf *.o $(OBJS) MagTFs
