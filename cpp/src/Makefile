srcfiles := main.o expData.o sysInfo.o
dstfiles := $(addprefix ../bin/,$(srcfiles))

a.out: main.o expData.o sysInfo.o energyExpressions.o twoParamsEnergies.o
	g++ ../bin/main.o ../bin/expData.o ../bin/sysInfo.o ../bin/energyExpressions.o ../bin/twoParamsEnergies.o  -o ../bin/a.out
	
main.o: main.cpp
	g++ -std=c++11 -o ../bin/main.o -c main.cpp

expData.o: expData.cpp ../include/expData.h
	g++ -std=c++11 -o ../bin/expData.o -c expData.cpp

sysInfo.o: sysInfo.cpp ../include/sysInfo.h
	g++ -std=c++11 -o ../bin/sysInfo.o -c sysInfo.cpp

energyExpressions.o: energyExpressions.cpp ../include/energyExpressions.h
	g++ -std=c++11 -o ../bin/energyExpressions.o -c energyExpressions.cpp

twoParamsEnergies.o: twoParamsEnergies.cpp ../include/twoParamsEnergies.h
	g++ -std=c++11 -o ../bin/twoParamsEnergies.o -c twoParamsEnergies.cpp
