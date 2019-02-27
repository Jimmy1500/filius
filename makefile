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
	${CXX} ${CXXFLAGS} -fPIC -o ${OBJ}$@ -c $<

G2PP.so: ${OBJ}G2PP.o ${OBJ}Simulation.o ${OBJ}Curve.o
	${CXX} ${CXXFLAGS} -o ${LIB}$@ -shared $^

Optimization.so: ${OBJ}Optimization.o ${OBJ}Swaption.o ${LIB}G2PP.so
	${CXX} ${CXXFLAGS} -o ${LIB}$@ -shared $^

TestSimulation: ${OBJ}TestSimulation.o ${OBJ}Swaption.o ${LIB}G2PP.so ${LIB}Optimization.so
	${CXX} ${CXXFLAGS} -o ${BIN}$@ $^

.PHONY: clean build run
clean:
	( rm -f ${OBJ}* ${LIB}* ${BIN}* )

build:
	( mkdir -p ${OBJ} ${LIB} ${BIN} && make )

run:
	${BIN}${EXECUTABLES}
