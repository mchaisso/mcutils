#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
using namespace std;

void PrintUsage() {
	cout << "Usage: samToBed file.sam [--reportAccuracy]" << endl;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		PrintUsage();
		exit(1);
	}

	ifstream samIn(argv[1]);
	int argi = 2;
	bool reportAccuracy = false;

	while (argi < argc) {
		if (strcmp(argv[argi], "--reportAccuracy") == 0) {
			reportAccuracy = true;
		}
		else {
			PrintUsage();
			cout << "Unknown option " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}

	while (samIn) {
		string line;
		getline(samIn, line);
		
		if (line.size() > 0 and line[0] == '@') {
			continue;
		}
		if (line.size() == 0) {
			break;
		}

		stringstream lineStrm(line);
		string contig;
		int flag;
		string chrom;
		int pos;
		int mapq;
		string cigar;
		string t1,t2;
		int alnLength;
		string seq;
		int strandFlag = 0x16;
		lineStrm >> contig >> flag >> chrom >> pos >> mapq >> cigar >> t1 >> t2 >> alnLength >> seq;
		int strand = 0;
		if (flag & strandFlag) {
			strand = 1;
		}
		if (chrom == "*") {
			continue;
		}
		string cigarPrefix = cigar.substr(0,100);
		string cigarSuffix = cigar.substr(max(0, (int)cigar.size()-100), 100);
		int s1 = 0;
		int nMatch = 0, nIns = 0, nDel = 0, nMisMatch = 0;
		int frontClip = 0, endClip = 0;
		int i = 0;
		bool alnStarted = false;
		while (i < cigar.size()) {
			int len = atoi(&cigar[i]);
			while (i < cigar.size() and cigar[i] >= '0' and cigar[i] <= '9') { i++;}
			char op = cigar[i];
			if (op == 'S' or op == 'H') {
				if (op == 'S' and alnStarted == false) {
					frontClip = len;
				}
				else if (op == 'S' and alnStarted == true) {
					endClip = len;
				}
			}
			else {
				alnStarted = true;
				if (op == 'M') {
					nMatch += len;
				}
				if (op == '=') {
					nMatch +=len;
				}
				else if (op == 'X') {
					nMisMatch+=len;
				}
				else if (op == 'I') {
					nIns += len;
				}
				else if (op == 'D') {
					nDel += len;
				}
			}
			i+=1;
		}
		int denom = nMatch + nMisMatch + nIns + nDel;
		float identity = 0;
		if (denom > 0) {
			identity = ((float) nMatch) / denom;
		}
		cout << chrom << "\t" << pos - 1 << "\t" << pos + alnLength - 1  << "\t" 
				 << contig << "\t" << strand << "\t" << frontClip << "\t" << seq.size() - endClip << "\t" 
				 << mapq << "\t";
		if (reportAccuracy) {
			cout << identity << "\t"  
					 << nMatch << "\t" << nMisMatch << "\t"
					 << nIns << "\t" << nDel;
		}
		cout << endl;
	}
}
	
