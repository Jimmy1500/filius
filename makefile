OBJECTS=$(shell cd src && ls *.cpp | sed -e 's/cpp/o/g')
LIBS=G2PP.so \
	 Optimization.so
EXECUTABLES=TestSimulation
COMPILER=compiler.mk

include $(COMPILER)

# build targets
all: ${OBJECTS} ${LIBS} ${EXECUTABLES}

# build rules | $@: left(%), $<: right(first), $^: right(all)
Curve.o: ${SRC}Curve.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

Simulation.o: ${SRC}Simulation.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

G2PP.o: ${SRC}G2PP.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

Swaption.o: ${SRC}Swaption.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

Optimization.o: ${SRC}Optimization.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

TestSimulation.o: ${SRC}TestSimulation.cpp
	${CXX} ${CXXFLAGS} -fPIC -c $< -o ${OBJ}$@

G2PP.so: ${OBJ}G2PP.o ${OBJ}Simulation.o ${OBJ}Curve.o
	${CXX} ${CXXFLAGS} -shared $^ -o ${LIB}$@

Optimization.so: ${OBJ}Optimization.o ${LIB}G2PP.so ${OBJ}Swaption.o
	${CXX} ${CXXFLAGS} -shared $^ -o ${LIB}$@

TestSimulation: ${OBJ}TestSimulation.o ${OBJ}Swaption.o ${LIB}G2PP.so ${LIB}Optimization.so
	${CXX} ${CXXFLAGS} $^ -o ${BIN}$@

.PHONY: clean build run
clean:
	( rm -f ${OBJ}* ${LIB}* ${BIN}* )

build:
	( mkdir -p ${OBJ} ${LIB} ${BIN} && make )

run:
	${BIN}${EXECUTABLES}
