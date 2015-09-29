#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
using namespace std;
int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "Usage: samToBed file.sam " << endl;
		exit(1);
	}

	ifstream samIn(argv[1]);
	int tPos = 0;
	int qPos = 0;
	int contigIndex = 0;
	string previousContig = "";
	string previousChrom = "";
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

		lineStrm >> contig >> flag >> chrom >> pos >> mapq >> cigar >> t1 >> t2 >> alnLength >> seq;
		if ((previousContig != "" and contig != previousContig) or( previousChrom != "" and previousChrom != chrom)) {
			break;
		}
		if (chrom == "*") {
			continue;
		}
		int s1 = 0;
		int nMatch = 0, nIns = 0, nDel = 0, nMisMatch = 0;
		int frontClip = 0, endClip = 0;
		int i = 0;
		bool alnStarted = false;
		int seqLen = seq.size();
		int strand = flag && 0x16;
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
				int qPlotPos = qPos;
				if (strand != 0) {
					qPlotPos = seqLen - qPos;
				}
				if (op == 'M') {
					cout << tPos << "\t" << qPlotPos << "\t" << len << "\t" << contigIndex << "\t" << strand << endl;
					nMatch += len;
					tPos += len;
					qPos += len;
				}
				if (op == '=') {
					nMatch +=len;
					cout << tPos << "\t" << qPlotPos << "\t" << len << "\t" << contigIndex << "\t" << strand << endl;
					tPos += len;
					qPos += len;
				}
				else if (op == 'X') {
					nMisMatch+=len;
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
		previousChrom = chrom;
		contigIndex ++;
	}
}
	
