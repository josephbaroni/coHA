# ----------------------------
# Project: SE3GNN
# Directory structure:
# SE3GNN/
#   math/   (all .hpp and .cpp)
#   tests/  (all test .cpp)
# ----------------------------

CXX := g++-13
CXXFLAGS := -std=c++17 -O2 -Imath -Wall -Wextra -fPIC

# math source files
MATH_SRC := $(wildcard math/*.cpp)
MATH_OBJ := $(MATH_SRC:.cpp=.o)

# test source files
TEST_SRC := $(wildcard tests/*.cpp)

# library
LIB := libmath.a

# default target
all: $(LIB) test_equivariance

# build static library
$(LIB): $(MATH_OBJ)
	ar rcs $@ $^

# compile math objects
math/%.o: math/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# test executable
test_equivariance: $(TEST_SRC) $(LIB)
	$(CXX) $(CXXFLAGS) $^ -o $@

# clean
.PHONY: clean
clean:
	rm -f math/*.o *.o $(LIB) test_equivariance
