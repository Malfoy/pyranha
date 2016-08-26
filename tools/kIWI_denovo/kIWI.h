#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

struct pile{
    string kmer;
    vector<string> context;
    vector <uint> indexReads;
    int strand;
};


//~ struct comparePile{
    //~ bool operator()(const pile& p1, const edge& seqR){
        //~ return seqL.sequence <seqR.sequence;
    //~ }
//~ };

int findStrand(string s);
string getPrefix(uint n, const string& sequence);
string getSuffix(uint n, const string& sequence);
char revCompChar(char c);
vector<char> alt(char c);
string revComp(const string& s);
uint hammingDist(const string& s1, const string& s2, uint lim);
vector<pile> findExactPiles(pile& p);
void indexOverlaps(ifstream& readsFile, unordered_map <string, pile>& indexKmersPlus,  unordered_map <string, pile>& indexKmersMinus, vector <string>& sequences);
void detectPileForAStrand(unordered_map <string, pile>& indexKmers, unordered_map <string, vector <pile>>& newIndexKmers);
void writeResult(unordered_map <string, vector <pile>>& plus, unordered_map <string, vector <pile>>& minus, const vector <string>& sequences, const string& outFileName);
//~ string getCanonical(const string& seq);

