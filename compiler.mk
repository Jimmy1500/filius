# machine related info
PLATFORM=$(shell uname -s)
ARCHITECTURE=$(shell uname -m)
EXPORT=$(PLATFORM)_$(ARCHITECTURE)

# source tree directories
INC=include/
SRC=src/
OBJ=${EXPORT}/obj/
LIB=${EXPORT}/lib/
BIN=${EXPORT}/bin/

# compiler options
C        := gcc
CXX      := g++

DEBUG    := -g
CXX_STD  := -std=c++17
CXX_FLG  := --warn-unused-variable -unknown-pragmas -Wall -Werror
CXX_THR  := -pthread
CXX_OPT  := -Ofast
CXX_INC  := -I$(INC)
CXX_MAC  := -D__DEBUG__#-D__REGEN__

CFLAGS   := $(DEBUG) -gnu99 -Wall -Werror -pthread I$(INC)
CXXFLAGS += $(DEBUG)
CXXFLAGS += $(CXX_STD)
CXXFLAGS += $(CXX_FLG)
CXXFLAGS += $(CXX_THR)
CXXFLAGS += $(CXX_OPT)
CXXFLAGS += $(CXX_INC)
CXXFLAGS += $(CXX_MAC)
