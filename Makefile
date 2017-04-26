CONFIG_FILE := Makefile.config
ifneq ($(wildcard $(CONFIG_FILE)),)
include $(CONFIG_FILE)
endif

CC := mpic++
CFLAGS := -std=c++11 -g -Wall -O3
SRC_DIR := src
OBJ_DIR := build
SRCS := $(shell find $(SRC_DIR) ! -name "main.cc" ! -name "*test.cc" -name "*.cc")
HEADERS := $(shell find $(SRC_DIR) -name "*.h")
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
MAIN := src/main.cc
EXE := hci

.PHONY: all clean

all: $(EXE)

clean:
	rm -rf $(OBJ_DIR)/*
	rm -f ./hci

$(EXE): $(OBJS) $(MAIN) $(HEADERS)
	$(CC) $(CFLAGS) $(MAIN) $(OBJS) -o hci -lboost_mpi -lboost_serialization
	
$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CC) $(CFLAGS) -c $< -o $@
	
