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

DBG      := -g
CXX_STD  := -std=c++17
CXX_FLG  := --warn-unused-variable -unknown-pragmas -Wall -Werror
CXX_THD  := -pthread
CXX_OPT  := -Ofast
CXX_INC  := -I$(INC)
CXX_MCR  := -D__DEBUG__ -D__REGEN__

CFLAGS   := $(DBG) -gnu99 -Wall -Werror -pthread I$(INC)
CXXFLAGS += $(DBG)
CXXFLAGS += $(CXX_STD)
CXXFLAGS += $(CXX_FLG)
CXXFLAGS += $(CXX_THD)
CXXFLAGS += $(CXX_OPT)
CXXFLAGS += $(CXX_INC)
CXXFLAGS += $(CXX_MCR)
