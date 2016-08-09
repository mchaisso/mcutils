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
int Q = 0;
int T = 1;
int dir = Q;

class Block {
public:
	int t,q,l;
	Block(int &tp, int &qp, int &lp) {
		t = tp; q = qp; l=lp;
	}
	int Pos() const {
		if (dir == 0) {
			return q;
		}
		else {
			return t;
		}
	}
	int Map() {
		if (dir == 0) {
			return t;
		}
		else {
			return q;
		}
	}
	int operator<(const int p) const {
		return Pos() < p;
	}
};


typedef vector<Block> Blocks;
typedef vector<Blocks> MultiBlocks;
typedef vector<int> Strands;
typedef map<string, MultiBlocks> PosMap;
typedef map<string, Strands> StrandMap;
typedef map<string, vector<int> > LengthMap;
typedef map<string, vector<string> > ChromMap;


bool SetPos(Blocks &blocks, int index, int pos, int &mapPos) {
	assert(blocks[index].Pos() <= pos);
	if (index == blocks.size() - 1 and pos > blocks[index].Pos() + blocks[index].l) {
		return false;
	}
		
	if (index < blocks.size() - 1) {
		assert(blocks[index+1].Pos() >= pos);
	}
	if (blocks[index].Pos() + blocks[index].l >= pos) {
		mapPos = blocks[index].Map() + pos - blocks[index].Pos()  ;
	}
	else {
		// 
		// This position is in a gap.
		mapPos = blocks[index].Map() + blocks[index].l;
	}
	return true;
}

bool SearchBlocks(Blocks &blocks, int pos, int &mapPos) {
	if (blocks.size() == 0) {
		return 0;
	}
	assert(blocks.size() > 0);
	Blocks::iterator lb;
	lb = lower_bound(blocks.begin(), blocks.end(), pos);
	int i = lb - blocks.begin();
	//
	// Find the index that contains the block, or none if no such block exists.
	//
	if (i == blocks.size()-1) {
		int last = blocks.size()-1;
		if (blocks[last].Pos() + blocks[last].l >= pos  and blocks[last].Pos() < pos) {
			if (blocks[last].Pos() > pos) {
				i = last-1;		
			}
			else {
				i = last;
			}
		}
		else {
			return false;
		}
	}
	else if (i == 0) {
		//
		// Either the position is too early, or right at the start
		//
		if (blocks[i].Pos() == pos) {
			mapPos = blocks[i].Map();
			return true;
		}
		else {
			return false;
		}
	} else {
		i = i - 1;
	}

	assert(blocks[i].Pos() <= pos);
	return SetPos(blocks, i, pos, mapPos);
}

bool SearchContig(PosMap &posMap, ChromMap &chromMap, StrandMap &strandMap, LengthMap &lengthMap,
				  string contig, int pos,  string &mapChrom, int &mapPos, int &mapStrand) {
	if (posMap.find(contig) == posMap.end()) {
		return false;
	}
	int i;
	for (i = 0; i < posMap[contig].size(); i++) {
		int searchPos = pos;

		if (dir == Q) {
			if (strandMap[contig][i] != 0) {
				searchPos = lengthMap[contig][i] - pos - 1;
			}
		}
		if (SearchBlocks(posMap[contig][i], searchPos, mapPos) == true) {
			mapChrom = chromMap[contig][i];
			mapStrand = strandMap[contig][i];
			if (dir == T) {
				if (strandMap[contig][i] != 0) {
					mapPos = lengthMap[contig][i] - mapPos - 1;
				}
			}
			return true;
		}

	}
	return false;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: samLiftover file.sam coordinates.bed out.bed [options]" << endl;
		cout << "  --dir 0|1 (0)  Map from coordinates on the query to the target." << endl;
		cout << "                 A value of 1 maps from the target to the query." << endl;
		exit(1);
	}

	ifstream samIn(argv[1]);
	ifstream bedIn(argv[2]);
	ofstream bedOut(argv[3]);
	string outName = argv[3];
	
	outName += ".bad";
	ofstream badOut(outName.c_str());
	int argi = 4;
	bool printBedLine = false;
	if (argc >= 5) {
		while (argi < argc) {
			if (strcmp(argv[argi], "--dir") == 0) {
				dir = atoi(argv[++argi]);
			}
			if (strcmp(argv[argi], "--bedline") == 0) {
			  printBedLine = true;
			}
			++argi;
		}
	}
		
	if (bedIn.good() == false) {
		cerr << "Could not open " << argv[2] << endl;
		exit(1);
	}
	if (samIn.good() == false) {
		cerr << "Could not open " << argv[1] << endl;
		exit(1);
	}		

	int contigIndex = 0;
	string previousContig = "";
	map<string, MultiBlocks > posMap;
	map<string, Strands> strands;
	LengthMap lengths;
	map<string, vector<string> > chromMap;
	cerr << "Building map database." << endl;

	while (samIn) {
		string line;
		getline(samIn, line);
		int tPos = 0;
		int qPos = 0;

		if (line.size() > 0 and line[0] == '@') {
			continue;
		}
		if (line.size() == 0) {
			break;
		}
		stringstream lineStrm(line);
		string contig;
		string flagStr;
		string chrom;
		int pos;
		int mapq;
		string cigar;
		string t1,t2;
		int alnLength;
		string seq;
		lineStrm >> contig >> flagStr >> chrom >> pos >> mapq >> cigar >> t1 >> t2 >> alnLength >> seq;
		int flag;
		flag = strtol(flagStr.c_str(), NULL, 16);

		if (chrom == "*") {
			continue;
		}
		int s1 = 0;
		int nMatch = 0, nIns = 0, nDel = 0, nMisMatch = 0;
		int frontClip = 0, endClip = 0;
		int i = 0;
		bool alnStarted = false;

		string fromChrom = contig;
		string toChrom   = chrom;
		if (dir == T) {
			fromChrom = chrom;
			toChrom   = contig;
		}
		tPos = pos-1;
		int slash = fromChrom.find('/');
		if (slash != fromChrom.npos) {
		  fromChrom = fromChrom.substr(0,slash);
		}

		if (posMap.find(fromChrom) == posMap.end()) {

			posMap[fromChrom] = MultiBlocks();
			lengths[fromChrom] = vector<int>();
			chromMap[fromChrom] = vector<string>();
		}
		lengths[fromChrom].push_back(seq.size());
		posMap[fromChrom].push_back(Blocks());
		strands[fromChrom].push_back(flag & 0x16);
		
		chromMap[fromChrom].push_back(toChrom);
		Blocks *curBlocks = &(posMap[fromChrom][posMap[fromChrom].size()-1]);

		while (i < cigar.size()) {
			int len = atoi(&cigar[i]);
			while (i < cigar.size() and cigar[i] >= '0' and cigar[i] <= '9') { i++;}
			char op = cigar[i];
			if (op == 'S' or op == 'H') {
				if (op == 'S' and alnStarted == false) {
					frontClip = len;
					qPos = frontClip;
				}
				else if (op == 'S' and alnStarted == true) {
					endClip = len;
				}
			}
			else {
				alnStarted = true;
				if (op == 'M' or op == '=' or op == 'X') {
					curBlocks->push_back(Block(tPos, qPos, len));
					nMatch += len;
					tPos += len;
					qPos += len;
				}
				else if (op == 'I') {
					nIns += len;
					qPos += len;
				}
				else if (op == 'D') {
					nDel += len;
					tPos += len;
				}
			}
			i+=1;
		}
		previousContig = contig;
		contigIndex ++;
	}
	
	string bedLine;

	while (getline(bedIn, bedLine)) {
		if (bedLine.size() == 0) {
			break;
		}
		stringstream bedStrm(bedLine);
		string chrom, mapChrom, mapEndChrom;
		int start, mapStart, end, mapEnd, mapFrontStrand, mapEndStrand;
		bedStrm >> chrom >> start >> end;
		vector<string> remainder;
		while (bedStrm) {
		  string val;
		  bedStrm >> val;
		  if (val != "") {
			remainder.push_back(val);
		  }
		}
		bool foundStart, foundEnd;
		
		foundStart = SearchContig(posMap, chromMap, strands, lengths,
								  chrom, start, mapChrom, mapStart, mapFrontStrand);
		foundEnd = SearchContig(posMap, chromMap, strands, lengths,
								chrom, end, mapEndChrom, mapEnd, mapEndStrand);

		if (mapFrontStrand == mapEndStrand and mapFrontStrand != 0) {
		  int temp = mapStart;
		  mapStart = mapEnd; 
		  mapEnd   = temp;
		}

		if (foundStart == false or foundEnd == false or mapChrom != mapEndChrom) {
		  badOut << bedLine << " " << (int) foundStart << " " << (int) foundEnd << " " << (int) (mapStart == mapEnd) << endl;
		}
		else {
			bedOut << mapChrom << "\t" << mapStart << "\t" << mapEnd;
			int s;
			for (s = 0; s < remainder.size(); s++) {
				bedOut << "\t" << remainder[s];
			}
			//
			// Add the originals just in case
			//
			bedOut << "\t" << chrom << "\t" << start << "\t"<< end;
			
			bedOut << endl;
		}
	}
	
}
	
