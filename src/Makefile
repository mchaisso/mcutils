TARGETS=samToBed samToDotPlot samLiftover
INSTALLDIR?=../bin
all: $(TARGETS)
install:
	mkdir -p $(INSTALLDIR)
	cp $(TARGETS) $(INSTALLDIR)

CPPOPTS= -g
samToBed: SamToBed.cpp
	g++ $(CPPOPTS) $^ -o $@

samToDotPlot: SamToDotPlot.cpp
	g++ $(CPPOPTS) $^ -o $@

samLiftover: SamLiftover.cpp
	g++ $(CPPOPTS) $^ -o $@


clean:
	/bin/rm $(TARGETS)


