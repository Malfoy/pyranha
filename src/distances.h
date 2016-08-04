#ifndef DISTANCE
#define DISTANCE
#include <string>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cctype>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "utils.h"
#include "ssw_cpp.h"


using namespace std;


uint distHamming(const string& read, const string& reference, uint maxMissmatch);
uint distHammingIndel(const string& read, const string& reference, uint maxMissmatch);
void printAlignmentSW(const StripedSmithWaterman::Alignment& alignment);
void alignSW(const string& ref, const string& query);
int32_t nbMismatchesSW(const string& ref, const string& query);



#endif
