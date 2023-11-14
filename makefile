CXX					:=  g++
CXXFLAGS				:=  -std=c++11 -O3 -fopenmp

INCLUDES				:= $(wildcard src/*.h)
SRCS					:= $(wildcard src/*.cpp) 
OBJS					:= $(patsubst %.cpp, %.o, $(SRCS))

MagTF_INCLUDE		:= -I src/
EIGEN_INCLUDE		:= -I contrib/eigen

%.o : %.cpp 
	@echo "MagTF is compiling "$<"..."
	@$(CXX) $(CXXFLAGS) $(MagTF_INCLUDE) $(EIGEN_INCLUDE) -c $< -o $@

MagTF: $(OBJS) $(INCLUDES)
	@$(CXX) $(CXXFLAGS) -o MagTF $(OBJS)
	
clean:	
	@rm -rf *.o $(OBJS) MagTF
