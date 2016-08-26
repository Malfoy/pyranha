#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include "kIWI.h"

using namespace std;


/* detect strand of the read */
int findStrand(string s){
	size_t plus= s.find("+");
	if (plus != string::npos){
		return 1;
	} else {
		size_t minus= s.find("-");
		if (minus != string::npos){
			return 0;
		} else {
			return -1;
		}
	}
}

/* reverse complement of a nt */
char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}

/* reverse complement of a sequence */
string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}


vector<char> alt(char c) {
	switch (c) {
		case 'A': return {'C','G','T'};
		case 'C': return {'A','G','T'};
		case 'G': return {'A','C','T'};
	}
	return {'A','C','G'};
}

//~ string getCanonical(const string& seq){
	//~ return min(seq, revComp(seq));
//~ }


/* n-prefix of a sequence */
string getPrefix(uint n, const string& sequence){
	return sequence.substr(0, n);
}


/* n-suffix of a sequence */
string getSuffix(uint n, const string& sequence){
	return sequence.substr(sequence.size()-n, n);
}



/* index 10-mers overlaps found in case of ping pong. We index in a map: 10mers -> Pile. A Pile struct contains the strands and read indexes*/
void indexOverlaps(ifstream& readsFile, unordered_map <string, pile>& indexKmersPlus,  unordered_map <string, pile>& indexKmersMinus, vector <string>& sequences){
	string line;
	uint indexRead(0);
	while (not readsFile.eof()){
        getline(readsFile, line);  // get header
		int strand(findStrand(line));
		if (strand != -1){  // we found a strand in the header
			if (not line.empty()){
				getline(readsFile, line);
				string overlap, overlap_temp, context;
				if (strand){  // plus strand
					overlap_temp = getPrefix(15, line);
					overlap = getPrefix(10, overlap_temp);
					context = getSuffix(5, overlap_temp);
					vector <string> vContext({context});
					vector <uint> vRead({indexRead});
					if (indexKmersPlus.count(overlap)){
						indexKmersPlus[overlap].context.push_back(context);
						indexKmersPlus[overlap].indexReads.push_back(indexRead);
					} else {
						pile pileP({overlap, vContext, vRead, strand});
						indexKmersPlus.insert({overlap, pileP});
					}
				} else {  // minus strand
					string rc(revComp(line));
					overlap_temp = getPrefix(15, rc);
					context = getSuffix(5, overlap_temp);
					overlap = getPrefix(10, overlap_temp);
					vector <string> vContext({context});
					vector <uint> vRead({indexRead});
					if (indexKmersMinus.count(overlap)){
						indexKmersMinus[overlap].context.push_back(context);
						indexKmersMinus[overlap].indexReads.push_back(indexRead);
					} else {
						pile pileP({overlap, vContext, vRead, strand});
						indexKmersMinus.insert({overlap, pileP});
					}
				}
				sequences.push_back(line);
				++ indexRead;
			}
		}
	}
}


/* compute Hamming distance between two strings */
uint hammingDist(const string& s1, const string& s2, uint lim){
	uint dist(0);
	if (s1 <= s2){
		for (uint i(0); i < s1.size(); ++i){
			if (s1[i] != s2[i]){
				++dist;
				if (dist >= lim){
					break;
				}
			}
		}
	} else {
		for (uint i(0); i < s2.size(); ++i){
			if (s1[i] != s2[i]){
				++dist;
				if (dist >= lim){
					break;
				}
			}
		}
	}
	return dist;
}


/* find piles allowing no sequencing errors/ variant */
vector<pile> findExactPiles(pile& p){
	vector<pile> singleton;
	vector<pile> newPiles;
	uint currentIndex(0), i(0);
	while (i < p.context.size()){
		vector <uint> vRead({p.indexReads[i]});
		vector <string> vContext({p.context[i]});
		pile toAdd({p.kmer, vContext, vRead, p.strand});
		if (i == 0){
			if (p.context[i] != p.context[i+1]){ // two different 5-mers
				singleton.push_back(toAdd);
			} else {  // similar 5-mers
				newPiles.push_back(toAdd);
			}
		} else if (i == p.context.size()-1){
			if (p.context[i] != p.context[i-1]){
				singleton.push_back(toAdd);
			} else {
				if (newPiles.empty()){
					newPiles.push_back(toAdd);
				} else { //  5-mer similar to previous : we only add the read index
					newPiles[currentIndex].indexReads.push_back(p.indexReads[i]);
				}
			}
		} else {
			if (p.context[i] != p.context[i+1] and p.context[i] != p.context[i-1]){
				singleton.push_back(toAdd);
				if (currentIndex != 0){
					++currentIndex; // a pile has finished with this singleton, no other similar 5-mer will be found
				}
			} else {
				if (newPiles.empty()){
					newPiles.push_back(toAdd);
				} else { //  5-mer similar to previous : we only add the read index 
					newPiles[currentIndex].indexReads.push_back(p.indexReads[i]);
				}
			}
		}
		++i;
	}
	newPiles.insert(newPiles.end(), singleton.begin(), singleton.end());
	return newPiles;
}



//~ void detectPileForAStrand(unordered_map <string, pile>& indexKmers, unordered_map <string, vector <pile>>& newIndexKmers){
	//~ cout << indexKmers.size() << endl;
	//~ for (auto i(indexKmers.begin()); i!= indexKmers.end(); ++i){
		//~ vector <pile> newPile;
		//~ uint pileIndex;
		//~ if (i->second.context.size() > 2){ // a pile can exist
			//~ for (uint seq1(0); seq1 < second.context.size() - 1; ++seq1){
				//~ for (uint seq2(1); seq2 < second.context.size(); ++seq2){
				//~ uint dist(hammingDist(i->second.context[seq1], i->second.context[seq2], 2));
				//~ }
			//~ }
		//~ } else {
			//~ newPile.push_back(i->second);
		//~ }
		//~ newIndexKmers.insert({i->first, newPile});
	//~ }
//~ }

/* detects several rna molecule mapping on the same location and strand = one pile */
void detectPileForAStrand(unordered_map <string, pile>& indexKmers, unordered_map <string, vector <pile>>& newIndexKmers){
	for (auto i(indexKmers.begin()); i!= indexKmers.end(); ++i){
		vector <pile> newPile;
		if (i->second.context.size() > 2){ // a pile can exist
			sort(i->second.context.begin(), i->second.context.end());
			newPile = findExactPiles(i->second);
		} else {
			newPile.push_back(i->second);
		}
		newIndexKmers.insert({i->first, newPile});
	}
}



void writeResult(unordered_map <string, vector <pile>>& plus, unordered_map <string, vector <pile>>& minus, const vector <string>& sequences,const string& outFileName){
	ofstream out(outFileName);
	for (auto i(plus.begin()); i != plus.end(); ++i){
		if (minus.count(i->first)){ // overlap shared between two strands
			vector <uint> readsPlus, readsMinus;
			for (uint pileIndex(0); pileIndex < i->second.size(); ++pileIndex){
				out << "ping-pong: "<< i->first << endl;
				out << "+: ";
				for (uint j(0);  j< i->second[pileIndex].indexReads.size(); ++j){
					out <<  i->second[pileIndex].indexReads[j] << ", ";
				}
				out << endl;
			}
			for (uint pileIndex(0); pileIndex < minus[i->first].size(); ++pileIndex){
				out << "-: ";
				for (uint j(0);  j< minus[i->first][pileIndex].indexReads.size(); ++j){
					out <<  minus[i->first][pileIndex].indexReads[j] << ", ";
				}
				out << endl;
			}
		} else { //degenerated only
			for (uint pileIndex(0); pileIndex < i->second.size(); ++pileIndex){
				if (i->second[pileIndex].indexReads.size() > 1){
					out << "degenerated: "<< i->first << endl;
					out << "+: " ;
					for (uint j(0);  j< i->second[pileIndex].indexReads.size(); ++j){
						out <<  i->second[pileIndex].indexReads[j] << ", ";
					}
					out << endl;
				}
			}
		}
	}
	for (auto i(minus.begin()); i != minus.end(); ++i){
		if (not plus.count(i->first)){ // reverse only
			for (uint pileIndex(0); pileIndex < i->second.size(); ++pileIndex){
				if (i->second[pileIndex].indexReads.size() > 1){
					out << "reverse: "<< i->first << endl;
					out << "-: " ;
					for (uint j(0);  j< i->second[pileIndex].indexReads.size(); ++j){
						out <<  i->second[pileIndex].indexReads[j] << ", ";
					}
					out << endl;
				}
			}
		}
	}
}
