#include "utils.h"
#include "mapping.h"
#include <fstream>
#include <chrono>

#include <getopt.h>

using namespace std;


int main(int argc, char ** argv){
        srand(time(NULL));
        if(argc==1){
                cout<<"Usage : "<<endl;
                cout<<"-u : read file (testread.fa)"<<endl;
                cout<<"-x : reference file (testref.fa)"<<endl;
                cout<<"-b : kmer bank to index  (kmercount)"<<endl;
                cout<<"-k : size of anchor (15)"<<endl;
                cout<<"-m : max missmatch (2)"<<endl;
				cout<<"-t : thread number (1)"<<endl;
				cout<<"-f : index 1/fraction kmer (1)"<<endl;
                cout<<"-c : output unaligned reads sequence instead of 'not_aligned'"<<endl;

                return 0;
        }
        string refFile("testref.fa"),seq,ref,readFile("testread.fa"), tempAlign,kmerCountFile("kmercount");
        uint k(15);
        uint maxMiss(5);
		uint coreNumber(1);
		uint fraction(1);
        bool notAlignedSequence(false);
        char c;
        while ((c = getopt (argc, argv, "u:k:x:m:c:t:f:b:")) != -1){
                switch(c){
                        case 'u':
                                readFile=optarg;
                                break;
                        case 'x':
                                refFile=optarg;
                                break;
                        case 'b':
                                kmerCountFile=optarg;
                                break;
                        case 'k':
                                k=stoi(optarg);
                        break;
                        case 'm':
                                maxMiss=stoi(optarg);
                        break;
						case 't':
                                coreNumber=stoi(optarg);
                        break;
						case 'f':
                                fraction=stoi(optarg);
                        break;
                        case 'c':
                                tempAlign = optarg;
                                if (tempAlign == "TRUE" or tempAlign == "True" or tempAlign == "true" or tempAlign == "T" or tempAlign == "t"){
                                    notAlignedSequence = true;
                                }
                        break;
                }
        }
        auto startChrono2=chrono::system_clock::now();
        extractRef(refFile,ref);
        cout<<"Filling index of "<<refFile<<endl;
        unordered_map<kmer,vector<position>> kmer2pos;
        // fillIndex(refFile, k, kmer2pos,fraction);
		fillMPHF(ref,fraction,k,coreNumber,kmerCountFile);
		auto end2=chrono::system_clock::now();auto waitedFor2=end2-startChrono2;
		cout<<"Indexing took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor2).count())<<" sec"<<endl;
        cout<<"Mapping "<<readFile<<endl;
        auto startChrono=chrono::system_clock::now();
        uint nbread(mapReadFile(readFile,k,kmer2pos, ref,maxMiss, notAlignedSequence,coreNumber));
        auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
        cout<<"Mapping took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
        cout<<"Throughout: "<< nbread/(1000*(chrono::duration_cast<chrono::seconds>(waitedFor).count()))<<"k read by second or "
        << (nbread*3600/(1000000*(chrono::duration_cast<chrono::seconds>(waitedFor).count())))<<"M by hour"<<endl;
        return 0;
}
