OBJS=$(shell cd src && ls *.cpp | sed -e 's/cpp/o/g')
LIBS=G2PP.so 			\
	 Optimization.so
EXECUTABLES=TestSimulation
COMPILER=compiler.mk

include $(COMPILER)

# build targets
all: ${OBJS} ${LIBS} ${EXECUTABLES}

# build rules | $@: left(%), $<: right(first), $^: right(all)
%.o: $(SRC)%.cpp
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
