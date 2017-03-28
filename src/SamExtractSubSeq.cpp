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
using namespace std;
typedef pair<int,int> Coordinate;


#include "LiftOver.h"


int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: samExtractSubSeq file.sam out.fasta [options]" << endl;
		cout << "  --region   chr:start-end" << endl;
		cout << "  --bed      bed file of coordinates." << endl;
		exit(1);
	}
	int dir = Q;
	ifstream samIn(argv[1]);
	ofstream fastaOut(argv[2]);
	
	int argi = 3;
	bool printBedLine = false;
	vector<string> regions;
	string bedFileName = "";
	while (argi < argc) {
		if (strcmp(argv[argi], "--region") == 0) {
			regions.push_back(argv[++argi]);
		}
		else if (strcmp(argv[argi], "--bed") == 0) {
			bedFileName = argv[++argi];
		}
		else {
			cout << "Bad command: " << argv[argi] << endl;
			exit(1);
		}
		++argi;	
	}
	Block::dir = T;
	dir=T;
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
	if (bedFileName != "") {
		ifstream bedIn;
		bedIn.open(bedFileName.c_str());
		
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
			
			fastaOut << ">" << chroms[r] << ":" << starts[r] << "-" << ends[r] << "/" << mapChrom << ":" << mapStart << "-" << mapEnd << endl;

			string seq = seqMap[startContig][startContigIndex].substr(seqStart, seqEnd - seqStart);
			fastaOut << seq << endl;

			
		}
	}
	
}
	
