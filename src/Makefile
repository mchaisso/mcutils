TARGETS=samToBed samToDotPlot samLiftover subdepth samSubseq
INSTALLDIR?=../bin
all: $(TARGETS)
install:
	mkdir -p $(INSTALLDIR)
	cp $(TARGETS) $(INSTALLDIR)

CPPOPTS= -g
subdepth: SubSampleDepth.cpp
	g++ $(CPPOPTS) $^ -o $@

samToBed: SamToBed.cpp
	g++ $(CPPOPTS) $^ -o $@

samToDotPlot: SamToDotPlot.cpp
	g++ $(CPPOPTS) $^ -o $@

samLiftover: SamLiftover.cpp LiftOver.h
	g++ $(CPPOPTS) $< -o $@

samSubseq: SamExtractSubSeq.cpp LiftOver.h
	g++ $(CPPOPTS) $< -o $@


clean:
	/bin/rm $(TARGETS)


