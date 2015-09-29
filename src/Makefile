TARGETS=samToBed samToDotPlot samLiftover

all: $(TARGETS)
install:
	cp $(TARGETS) ../bin

CPPOPTS= -O3
samToBed: SamToBed.cpp
	g++ $(CPPOPTS) $^ -o $@

samToDotPlot: SamToDotPlot.cpp
	g++ $(CPPOPTS) $^ -o $@

samLiftover: SamLiftover.cpp
	g++ $(CPPOPTS) $^ -o $@


clean:
	/bin/rm $(TARGETS)


