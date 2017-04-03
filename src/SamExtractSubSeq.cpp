#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <map>
#include <algorithm>
#include <vector>
#include <assert.h>

#include <stdio.h>
#include <unistd.h>


using namespace std;
typedef pair<int,int> Coordinate;


#include "LiftOver.h"

struct opts{
    std::string outputfile;
    std::string samfile   ;
    std::string prefix    ;
    std::string region    ;
    std::string bedfile   ;
}globalOptions;


void printHelp(void){
    std::stringstream helpMessage;
    helpMessage << "Usage: samExtractSubSeq -s file.sam -f out.fasta -r chr:start-end" << std::endl;
    helpMessage << "       -r - optional - <STRING> A target position to extract" << std::endl;
    helpMessage << "       -b - optional - <STRING> A bedfile of regions " << std::endl;
    helpMessage << "       -s - required - <STRING> The SAM file containing query to target alignments" << std::endl;
    helpMessage << "       -f - optional - <STRING> The output fasta [STDOUT]" << std::endl;
    helpMessage << "       -p - optional - <STRING> A prefix for the output sequence name [NONE]" << std::endl;
    std::cerr << helpMessage.str() << std::endl;
};

int main(int argc, char* argv[]) {

    int c;
    while ((c = getopt(argc, argv, "hr:s:f:p:b:")) != -1){
        switch (c){
        case 'h':
            {
                printHelp();
                exit(1);
            }
        case 'r':
            {
                globalOptions.region = optarg;
                std::cerr << "INFO: setting region: " << globalOptions.region << std::endl;
                break;
            }
        case 's':
            {
                globalOptions.samfile = optarg;
                break;
            }
        case 'f':
            {
                globalOptions.outputfile = optarg;
                break;
            }
        case 'p':
            {
                globalOptions.prefix = optarg;
                break;
            }
        case 'b':
            {
                globalOptions.bedfile = optarg;
            }
        case '?':
            {
            break;
            }
        }
    }

    if(globalOptions.samfile.empty()){
        printHelp();
        exit(1);
    }
	ifstream samIn(globalOptions.samfile.c_str());

    if (samIn.good() == false) {
        cerr << "FATAL: could not open sam file: " << argv[1] << endl;
        exit(1);
    }


    if(globalOptions.region.empty() && globalOptions.bedfile.empty()){
        std::cerr << "FATAL: need to specify a region -r OR -b a bedfile" << std::endl << std::endl;
        exit(1);
    }


	int dir = Q;
    Block::dir = T;
    dir=T;

    ofstream fastaOut;
    if(! globalOptions.outputfile.empty()){
        fastaOut.open(globalOptions.outputfile.c_str());
    }

	bool printBedLine = false;
	vector<string> regions;

    if(!globalOptions.region.empty()){
        regions.push_back(globalOptions.region);
    }

	if (samIn.good() == false) {
		cerr << "Could not open " << argv[1] << endl;
		exit(1);
	}

	int contigIndex = 0;
	map<string, MultiBlocks > posMap;
	map<string, Strands> strands;
	LengthMap lengths;
	map<string, vector<string> > chromMap;
	map<string, vector<string> > seqMap;
	typedef map<string, vector<int>  > ClipMap;
	ClipMap clipMap;
	cerr << "Building map database." << endl;
	int nContigs = 0;
    nContigs = BuildMapDB(samIn, dir, posMap, strands, lengths, chromMap, seqMap, clipMap, true);
	map<string, Strands>::iterator strandIt, strandEnd;
	for (strandIt = strands.begin(); strandIt != strands.end(); ++strandIt) {
		int si;
		for (si = 0; si < (*strandIt).second.size(); si++) {
			(*strandIt).second[si] = 0;
		}
	}

	string bedLine;

	int r;
	vector<string> chroms;
	vector<int> starts;
	vector<int> ends;
	for (r = 0; r < regions.size(); r++) {

    size_t start_pos = 0;
    while((start_pos = regions[r].find(",", start_pos)) != std::string::npos) {
        regions[r].replace(start_pos, 1, "");
        start_pos += 0; // Handles case where 'to' is a substring of 'from'
    }
		size_t chrPos = regions[r].find(":", 0);
		if (chrPos == std::string::npos) {
			cout << "Malformatted region " << regions[r];
			exit(1);
		}
		else {
			chroms.push_back(regions[r].substr(0,chrPos));
		}
		size_t startPos=chrPos+1;
		size_t hyphen=regions[r].find("-", startPos);
		if (hyphen == std::string::npos) {
			cout << "Malformatted region " << regions[r];
			exit(1);
		}
		else {
			starts.push_back(atoi(regions[r].substr(startPos, hyphen - startPos).c_str()));
			ends.push_back(atoi(regions[r].substr(hyphen+1).c_str()));
		}
	}
	if (!globalOptions.bedfile.empty()) {
		ifstream bedIn;
		bedIn.open(globalOptions.bedfile.c_str());

		while (getline(bedIn, bedLine)) {
			if (bedLine.size() == 0) {
				break;
			}
			else {
				stringstream strm(bedLine);
				string chrom;
				int s;
				int e;
				strm >> chrom >> s >> e;
				chroms.push_back(chrom);
				starts.push_back(s);
				ends.push_back(e);
			}
		}
	}


	for (r = 0; r < chroms.size(); r++ ) {

		string chrom, mapChrom, mapEndChrom;
		int start, mapStart, end, mapEnd, mapFrontStrand, mapEndStrand;

		vector<string> remainder;

		bool foundStart=false, foundEnd=false;
		string startContig="", endContig="";
		int startContigIndex=0, endContigIndex=0;
		foundStart = SearchContig(posMap, chromMap, strands, lengths,
															chroms[r], starts[r], mapChrom, mapStart, mapFrontStrand, startContig, startContigIndex);


        // end-1 is the 0-based coordinate for the end of the alignment
        // and SearchContig uses 0-based coordinates (DG, 161205)
		if (foundStart) {
			//
			// Do search for next position, but starting in the block of the previous position.
			//
			foundEnd = SearchContig(posMap, chromMap, strands, lengths,
															chroms[r], ends[r], mapEndChrom, mapEnd, mapEndStrand, endContig, endContigIndex, startContigIndex);
		}

		if (foundStart == true and
				foundEnd == true and
				mapChrom == mapEndChrom and
				startContig == endContig ) {
			if (seqMap.find(startContig) == seqMap.end() or
					seqMap[startContig].size() <= startContigIndex or
					startContigIndex != endContigIndex ) {
				cerr << "WARNING: Looup failed for " << chroms[r] << ":" << starts[r] << "-" << ends[r] << "/" << mapChrom << ":" << mapStart << "-" << mapEnd << endl;
				continue;
			}

			int clipStart= 0;

			if (clipMap.find(startContig) != clipMap.end() and clipMap[startContig].size() > startContigIndex) {
				clipStart = clipMap[startContig][startContigIndex];
			}
			int seqStart = mapStart - clipStart;
			int clipEnd=0;
			if (clipMap.find(startContig) != clipMap.end() and clipMap[startContig].size() > startContigIndex) {
				clipEnd = clipMap[startContig][startContigIndex];
			}

			int seqEnd  = mapEnd   - clipEnd;

			assert(seqStart <= seqEnd);
			if (seqStart == seqEnd) {
				continue;
			}
			if (seqStart >= seqMap[startContig][startContigIndex].size() or
					seqEnd  >= seqMap[startContig][endContigIndex].size()) {
				cerr << "WARNING: Looup failed for " << chroms[r] << ":" << starts[r] << "-" << ends[r] << "/" << mapChrom << ":" << mapStart << "-" << mapEnd << endl;
				continue;
			}

            std::stringstream ss;

			ss << ">" << globalOptions.prefix << chroms[r] << ":" << starts[r] << "-" << ends[r] << "/" << mapChrom << ":" << mapStart << "-" << mapEnd << endl;
			string seq = seqMap[startContig][startContigIndex].substr(seqStart, seqEnd - seqStart);

            if(!globalOptions.outputfile.empty()){
                fastaOut << ss.str();
                fastaOut << seq << std::endl;
            }
            else{
                std::cout << ss.str();
                std::cout << seq << std::endl;
            }
		}
	}
	return 0;
}

