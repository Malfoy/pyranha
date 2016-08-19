#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "utils.h"
#include "ssw_cpp.h"


using namespace std;

ofstream out("debug.txt");

uint distHamming(const string& read, const string& reference, uint maxMissmatch){
    uint error(0);
    for(uint i(0);i<read.size();++i){
        if(reference[i]==':'){return maxMissmatch+1;}
        if(read[i]!=reference[i]){
            if(++error>maxMissmatch){
                return error;
            }
        }
    }
    return error;
}


uint distHammingIndel(const string& read, const string& reference, uint maxMissmatch){
    uint error(0);
    for(uint i(0);i<read.size();++i){
        if(reference[i]==':'){return maxMissmatch+1;}
        if(read[i]!=reference[i]){
            if(i==read.size()-1){return ++error;}
            if(read[i+1]==reference[i]){//insertion
                if(i+2>=read.size() or i+1>=reference.size() or ++error>maxMissmatch){
                    return error;
                }
                return error + distHamming(read.substr(i+2), reference.substr(i+1), maxMissmatch-error);
            }else if(read[i]==reference[i+1]){//deletion
                if(i+2>=reference.size() or i+1>=read.size() or ++error>maxMissmatch){
                    return error;
                }
                return error + distHamming(read.substr(i+1), reference.substr(i+2), maxMissmatch-error);
            } else {//missmatch
                if(++error>maxMissmatch){
                    return error;
                }
            }
        }
    }
    // out<<"ret"<<endl;;
    return error;
}


void printAlignmentSW(const StripedSmithWaterman::Alignment& alignment){
  cout << "===== SSW result =====" << endl;
  cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
       << "Reference start:\t" << alignment.ref_begin << endl
       << "Reference end:\t" << alignment.ref_end << endl
       << "Query start:\t" << alignment.query_begin << endl
       << "Query end:\t" << alignment.query_end << endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
       << "Number of mismatches:\t" << alignment.mismatches << endl
       << "Cigar: " << alignment.cigar_string << endl;
  cout << "======================" << endl;
}



void alignSW(const string& ref, const string& query){
  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);
  printAlignmentSW(alignment);
}


int32_t nbMismatchesSW(const string& ref, const string& query){
  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);
  return alignment.mismatches;
}
