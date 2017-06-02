# Basic configurations.
CXX := g++
CXXFLAGS := -std=c++14 -Wall -Wextra -O3
SRC_DIR := src
BUILD_DIR := build
OUT := hci.out

# Intermediate variables.
MAIN := $(SRC_DIR)/main.cc
SRCS := $(shell find $(SRC_DIR) ! -name "main.cc" ! -name "*_test.cc" -name "*.cc")
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(BUILD_DIR)/%.o)
HEADERS := $(shell find $(SRC_DIR) -name "*.h")

.PHONY: all clean

all: $(OUT)

clean:
	rm -rf $(BUILD_DIR)
	rm $(OUT)

$(OUT): $(OBJS) $(MAIN) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(MAIN) -o $(OUT) $(LDLIBS)

$(OBJS): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CXX) $(CXXFLAGS) -c $< -o $@