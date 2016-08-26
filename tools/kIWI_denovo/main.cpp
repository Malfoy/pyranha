#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "kIWI.h"

using namespace std;

int main(int argc, char ** argv){
	string readFileName = argv[1];
	string outFile = argv[2];
	ifstream readFile(readFileName);
	unordered_map <string, pile> indexKmersPlus;
	unordered_map <string, pile> indexKmersMinus;
	vector <string> sequences;
	cout << "Indexing overlaps..." << endl;
	indexOverlaps(readFile, indexKmersPlus , indexKmersMinus, sequences);
	cout << "Detecting piles..." << endl;
	unordered_map <string, vector <pile>> newIndexKmersPlus;
	unordered_map <string, vector <pile>> newIndexKmersMinus;
	detectPileForAStrand(indexKmersPlus, newIndexKmersPlus);
	detectPileForAStrand(indexKmersMinus, newIndexKmersMinus);
	cout << "Writing output..." << endl;
	writeResult(newIndexKmersPlus, newIndexKmersMinus, sequences, outFile);
	return 0;
}
