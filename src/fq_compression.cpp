// vim: noai:ts=2:sw=2
//
//  FASTQcompression.cpp
//  FASTQcompression
//
//  Created by Giovanna on 28/05/20.
//  Copyright © 2020 Giovanna. All rights reserved.
//

#include <iterator>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include <stack>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include "include.hpp"
#include "dna_string_n.hpp"
#include "dna_bwt_n.hpp"
#include "../external/rankbv/rankbv.hpp"

#ifndef DEBUG
  #define DEBUG 0
#endif

#ifndef REVC 
  #define REVC 0
#endif

#define LONGEST 10000 //longest read

/*uint64_t min_cluster_length = 1000;
uint64_t num_close_freq = 0;
uint64_t num_cluster_5 = 0;
uint64_t cluster_5_close_freq = 0;
uint64_t num_cluster = 0;
uint64_t tryz = 0;*/

using namespace std;

string input_dna;
string input_qual;
string original_fastq;
string output;

bool ignore_headers = false; //ignore headers

/*
 * Debug mode variables: print the BWT and read names/qualities for each base
 */

bool debug = false; //print debug info
bool verbose = false; //verbose output
int max_id_len=20; //in read_info, store at most this number of chars for the IDs
vector<string> read_info;//if debug, store read coordinate for every BWT position
//------------

vector<string> read_ids;//ID of each read

uint64_t modified = 0;//count how many bases have been modified
uint64_t clusters_size=0;//total number of bases inside clusters

//vector<uint64_t> freqs(256,0);//temporary vector used to count frequency of bases inside the currently analyzed cluster
//vector<char> high_freqs;

vector<uint64_t> statistics_qual_before(256,0);//count absolute frequencies of qualities in reads, before modifying
vector<uint64_t> statistics_qual_after(256,0);//count absolute frequencies of qualities in reads, after modifying

int border = 1;//exclude/include this number of bases at the borders of the analyzed cluster

//minimum LCP required inside clusters
int K_def = 16;
int K = 0;

//do not consider clusters smaller than this threshold
int m_def = 1;
int m = 0;

#if REVC
  //this flag records if the input read file is divided in two halves: reads and their reverse complement
  bool revc = false;
#endif

//terminator character at the end of the reads
char TERM = '#';

string QUAL;//string of length |BWT| that contains the base qualities, for each BWT position
//string BWT_MOD;//string of length |BWT| that duplicates the BWT.
//char *BWT_MOD;
string BWT_MOD;
rankbv_t* rbv = NULL;

vector<bool> LCP_minima;//bitvector that stores the LCP minima
vector<bool> LCP_threshold;//bitvector that stores LCP values that exceed the threshold: >= K
//bool *LCP_minima;//bitvector that stores the LCP minima
//bool *LCP_threshold;//bitvector that stores LCP values that exceed the threshold: >= K

dna_bwt_n_t bwt;//the BWT data structure

float rare_threshold = 40;//Thresholds used to determinate which bases to discard from the cluster
float quality_threshold = 20; 

char default_value = '5';

char DNA[5] = {'A','C','G','T','N'};
#define dna(i) (DNA[i])
              //A B C D E F G H I J K L M N O P Q R S T U V W X Y
int ORD[25] = {0,0,1,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,3,0,0,0,0,0};
#define ord(c) (ORD[c-65])

void help(){
    
    cout << "FASTQcompression [options]" << endl <<
    "Options:" << endl <<
    "-h          Print this help." << endl <<
    "-e <arg>    Input eBWT file (A,C,G,T,#) of DNA (REQUIRED)." << endl <<
    "-q <arg>    Qualities permuted according to the DNA's ebwt (REQUIRED)." << endl <<
    "-f <arg>    Original fastq file (REQUIRED)." << endl <<
    "-o <arg>    Output fastq (REQUIRED)." << endl <<
    "-r          The second half of the reads is the reverse complement of the first half." << endl <<
    "-k <arg>    Minimum LCP required in clusters. Default: " << K_def << "." << endl <<
    "-m <arg>    Minimum length of cluster to be processed. Default: " << m_def << "." << endl <<
    "-t <arg>    ASCII value of terminator character. Default: " << int('#') << " (#)." << endl <<
    "-D          Print debug info for each BWT position." << endl << endl <<
    "-H          Ignores the headers." << endl <<
    
    "\nTo run FASTQcompression, you must first build the extended Burrows-Wheeler Transform " <<
    "of the input DNA sequences and the corresponding permutation of base quality scores." << endl;
    
    exit(0);
}

bool file_exists(string fileName){

  std::ifstream infile(fileName);

return infile.good();
}

void def_qs(){

  int min = 126;
  int max = 33;

	for(char& c : QUAL) {
		if((int)c != 0){
    	if((int)c > max){
				max = (int)c;
			}
			if((int)c < min){
				min = (int)c;
			}
		}
	}

  default_value = (char)(max+min)/2;

}


/*
 *  START PROCEDURES TO NAVIGATE SUFFIX TREE AND COMPUTE LCP MINIMA
 */

void update_LCP_leaf(sa_leaf L, uint64_t & lcp_values){
    
  for(uint64_t i = L.rn.first+1; i<L.rn.second; ++i){ 
    LCP_threshold[i] = (L.depth >= K);
    lcp_values++;
  }
}

void update_lcp_minima(sa_node_n x, uint64_t & n_min){
    
  /*
   * we have a minimum after the end of each child (that is different than #) of size at least 2 of the input node x, except
   * if the candidate minimum position is the last or exceeds the interval of x
   */
  
  if( x.first_C - x.first_A >= 2 and     // there are at least 2 'A'
     x.first_C < x.last-1             // candidate min in x.first_C is not >= last position
     ){
      LCP_minima[x.first_C] = true;
      n_min++;
  }
  
  if( x.first_G - x.first_C >= 2 and     // there are at least 2 'C'
     x.first_G < x.last-1             // candidate min in x.first_G is not >= last position
     ){
      LCP_minima[x.first_G] = true;
      n_min++;
  }

  if( x.first_N - x.first_G >= 2 and     // there are at least 2 'G'
     x.first_N < x.last-1             // candidate min in x.first_N is not >= last position
     ){
      LCP_minima[x.first_N] = true;
      n_min++;
  }

  if( x.first_T - x.first_N >= 2 and     // there are at least 2 'N'
     x.first_T < x.last-1             // candidate min in x.first_T is not >= last position
     ){
      LCP_minima[x.first_T] = true;
      n_min++;
  }
}

void detect_minima(){
    
  uint64_t n = bwt.size();
  
  cout << "\nPhase 2/4: navigating suffix tree leaves." << endl;
  
  /*
   * LCP_threshold[i] == 1 iff LCP[i] >= K
   */
  LCP_threshold = vector<bool>(n,false);
  //LCP_threshold = new bool[n]{false};
  
  uint64_t leaves = 0;//number of visited leaves
  uint64_t max_stack = 0;
  uint64_t lcp_values = 1;//number of computed LCP values
  
  {
    auto TMP_LEAVES = vector<sa_leaf>(5);
    
    stack<sa_leaf> S;
    S.push(bwt.first_leaf());
    
    int last_perc_lcp = -1;
    int perc_lcp = 0;
    
    while(not S.empty()){
        
      sa_leaf L = S.top();
      S.pop();
      leaves++;
      
      assert(leaf_size(L)>0);
      max_stack = S.size() > max_stack ? S.size() : max_stack;
      
      update_LCP_leaf(L,lcp_values);
      
      int t = 0;//number of children leaves
      bwt.next_leaves(L, TMP_LEAVES, t, 2);
      
      for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);
      perc_lcp = (100*lcp_values)/n;
      
      if(perc_lcp > last_perc_lcp){
        #if DEBUG
          if(verbose) cout << "LCP: " << perc_lcp << "%." <<  endl;
        #endif  
        last_perc_lcp = perc_lcp;
      }
    }
  }
  
  cout << "Computed " << lcp_values << "/" << n << " LCP threshold values." << endl;
  
  cout << "Max stack depth = " << max_stack << endl;
  cout << "Processed " << leaves << " suffix-tree leaves." << endl << endl;
  
  cout << "Phase 3/4: computing LCP minima." << endl;
  
  LCP_minima = vector<bool>(n,false);
  //LCP_minima = new bool[n]{false};
  
  auto TMP_NODES = vector<sa_node_n>(5);
  
  uint64_t nodes = 0;//visited ST nodes
  max_stack = 0;
  
  stack<sa_node_n> S;
  S.push(bwt.root());
  
  int last_perc_lcp = -1;
  int perc_lcp = 0;
  uint64_t n_min = 0;//number of LCP minima
  
  while(not S.empty()){
      
      max_stack = S.size() > max_stack ? S.size() : max_stack;
      
      sa_node_n N = S.top();
      S.pop();
      nodes++;
      
      //compute LCP values at the borders of N's children
      update_lcp_threshold(N, LCP_threshold, lcp_values, K);
      
      update_lcp_minima(N, n_min);
      
      //follow Weiner links
      int t = 0;
      bwt.next_nodes(N, TMP_NODES, t);
      
      for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);
      
      perc_lcp = (100*lcp_values)/n;
      
      if(perc_lcp > last_perc_lcp){
          
          #if DEBUG
            if(verbose) cout << "LCP: " << perc_lcp << "%." << endl;
          #endif
          
          last_perc_lcp = perc_lcp;
          
      }
      
  }
  
  cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
  cout << "Found " << n_min << " LCP minima." << endl;
  cout << "Max stack depth = " << max_stack << endl;
  cout << "Processed " << nodes << " suffix-tree nodes." << endl;
  
  
}

/*
 *  END PROCEDURES TO NAVIGATE SUFFIX TREE AND COMPUTE LCP MINIMA
 */

//This function modifies a Quality Score to be Illumina 8 Level Binning compliant
int illumina_8_level_binning(int newqs){

  if (newqs >= 40) newqs = 40;
  else if (newqs >= 35) newqs = 37;
  else if (newqs >= 30) newqs = 33;
  else if (newqs >= 25) newqs = 27;
  else if (newqs >= 20) newqs = 22;
  else if (newqs >= 10) newqs = 15;
  else if (newqs >= 2) newqs = 6;

  return newqs+33;

}


//This function calculates the average quality score in a cluster
int avg_qs(uint64_t start, uint64_t end){

  int sum=0;
  int num=0;
  
  for(uint64_t j=start; j<=end; j++){
  
  	if(bwt[j] != TERM){
  		sum=sum+(int)QUAL[j];
  		num++;
  	}
  }

  if(sum==0) return 0;
  return (sum/num);
}


//This function calculates the max quality score in a cluster
int max_qs(uint64_t start, uint64_t end){

  int max=0;
  for(uint64_t j=start; j<=end; j++){
    if(bwt[j] != TERM){
      if((int)QUAL[j] > max){
        max = (int)QUAL[j];
      }
    }
  }
  return max;
}


//This function calculates the mean error and, relying on it, calculates a new quality score
int mean_error(uint64_t start, uint64_t end){

double avg_err=0;
double num = 0;
double sum_err = 0;

  for(uint64_t j=start; j<=end; j++){
    if(bwt[j] != TERM){
      num++;
      sum_err = sum_err + pow(10, -((double)QUAL[j]-33)/10);
    }
  }
  avg_err = sum_err/num;
  int qs = round(-10*log10(avg_err));

return qs+33;
}

/*
 *
 * START PROCEDURE TO ANALYZE A CLUSTER [begin, i]
 *
 * We may include/exclude some symbols at cluster beginning by changing the variable border
 */

void process_cluster(uint64_t begin, uint64_t i){

  uint64_t start=(begin>=border?begin-border:0);
  uint64_t end=(i>border?i-border:0);
  
  uint64_t size = (end-start+1);
  clusters_size += size;

  //cluster is too short
  if(size < m) return;

  char newqs;
  uint64_t maxfreq = 0;
  char mostfreq;
  uint64_t freqs[5]{0};
 
  
  //printing bases+QS in the cluster to look them up
  
  #if DEBUG
    if(verbose) cout << "----\n";
  #endif

  int base_num = 0;
  for(uint64_t j = start; j <= end; ++j){
    /*Counts the frequency of each base and stores it in a vector, moreover stores the maximum/avg QS in a variable*/
    if(bwt[j] != TERM){
      freqs[ord(bwt[j])]++;
      base_num++;
    }
    #if DEBUG
      if(verbose) cout << bwt[j] << "\t" << (int)QUAL[j]-33 << endl;
    #endif
  }
  
  if(base_num == 0) return;

  /*Through max_qs we obtain the highest qs in the cluster, through avg_qs we obtain the average qs in the cluster, while
  through default_value we set the new quality score value to a fixed value previously calculated */

  #if M==0
    newqs = max_qs(start,end);
  #elif M==1
    newqs = avg_qs(start,end);
  #elif M==2
    newqs = default_value;
  #elif M==3
    newqs = mean_error(start,end);
  #else
    cout << "WARNING: unsupported choice. The process will use M=0." << endl;
    newqs = max_qs(start,end);
  #endif
  #if B==1
    newqs = illumina_8_level_binning(newqs-33);
  #endif

  #if DEBUG
    if(verbose) cout << "****\n";
  #endif

  //Frequency of the bases in percentage	
  freqs[0] = (100*freqs[0])/(base_num); //'A'
  freqs[1] = (100*freqs[1])/(base_num); //'C'
  freqs[2] = (100*freqs[2])/(base_num); //'G'
  freqs[3] = (100*freqs[3])/(base_num); //'T'
  freqs[4] = (100*freqs[4])/(base_num); //'N'

  /*Through these variables we obtain the most frequent base in the cluster and its frequency */
  //mostfreq = std::max_element(freqs.begin(),freqs.end()) - freqs.begin();
  //maxfreq = *std::max_element(freqs.begin(), freqs.end());

  mostfreq = 0; //'A'
  if(freqs[1]>freqs[mostfreq]) mostfreq = 1;
  if(freqs[2]>freqs[mostfreq]) mostfreq = 2;
  if(freqs[3]>freqs[mostfreq]) mostfreq = 3;
  if(freqs[4]>freqs[mostfreq]) mostfreq = 4;

  maxfreq = freqs[mostfreq]; //freqs['A']
  mostfreq = dna(mostfreq);

  //Check if there's a base which frequence is similar to maxfreq and store it in a vector 
  /*
  vector<char> high_freqs;
  if(freqs['A'] >= maxfreq -10) high_freqs.push_back('A');
  if(freqs['C'] >= maxfreq -10) high_freqs.push_back('C'); 
  if(freqs['G'] >= maxfreq -10) high_freqs.push_back('G');
  if(freqs['T'] >= maxfreq -10) high_freqs.push_back('T');
  if(freqs['N'] >= maxfreq -10) high_freqs.push_back('N');
  */
  bool high_freqs[5]{false};
  int high_freqs_size=0;
  for(int i=0; i<5; i++){ //'A', 'C', 'G', 'T', 'N'
    if(freqs[i] >= maxfreq -10){
      high_freqs[i] = true; //high_freqs.push_back('A');
      high_freqs_size++;
    }
  }


  //If there are less than 5 bases and there's more than a base with high frequency we haven't enough information
  //to do corrections
  //if(high_freqs.size() > 1 && base_num < 5){
  if(high_freqs_size > 1 && base_num < 5){
    /*
    freqs['A'] = 0;
    freqs['C'] = 0;
    freqs['G'] = 0;
    freqs['T'] = 0;
    freqs['N'] = 0;
    high_freqs.clear();
    for(int i=0; i<5; i++){ //'A', 'C', 'G', 'T', 'N'
      freqs[i]=0;
      high_freqs[i]=false;
    }
    */
    return;
  }

  /*
  * If there's only one base with high frequency, in this cycle we modify the values of QS and, 
  * if the base is less frequent than rare_threshold and its QS is minor then quality_threshold, also the value stored in BWT_MOD
  */
  //if(high_freqs.size() == 1){
  if(high_freqs_size == 1){
    for(uint64_t j = start; j <= end; ++j){
      if(bwt[j] != TERM){
        if((freqs[ord(bwt[j])]) < rare_threshold){	
          if((int)(QUAL[j]-33) < quality_threshold){
            //BWT_MOD[j] = mostfreq;
            //FELIPE
            BWT_MOD.push_back(mostfreq);
            rankbv_setbit(rbv,j);
            modified++;
          }
        }
        //FELIPE
        /*
        if(newqs != QUAL[j]) cout<<QUAL[j]<<" --> "<<newqs<<endl;
        else cout<<QUAL[j]<<endl;
        */
        QUAL[j] = newqs;
      }
    }
  }
  //Otherwise we have to decide between more bases. 
  else{

    //To do so we use the backward search to find which character in EBWT precedes each base in the cluster.
    //vector<vector<uint64_t>> LF_vect(256, vector<uint64_t>(256,0));
    uint64_t LF_vect2[255][255]{0};

    for(uint64_t j = start; j <= end; ++j){
      if(bwt[j] != TERM){
        //LF_vect[bwt[j]][bwt[bwt.LF(j)]]++;
        LF_vect2[bwt[j]][bwt[bwt.LF(j)]]++; //bwt[bwt.LF(j)] can be '#'
      }
    }

    //pre-compute max's
    char LF_max[5]{0};
    for(int i=0; i<5; i++){
      int max=0;
      if(LF_vect2[dna(i)]['#']>LF_vect2[dna(i)][max]) max = '#'; 
      if(LF_vect2[dna(i)]['A']>LF_vect2[dna(i)][max]) max = 'A';
      if(LF_vect2[dna(i)]['C']>LF_vect2[dna(i)][max]) max = 'C';
      if(LF_vect2[dna(i)]['G']>LF_vect2[dna(i)][max]) max = 'G';
      if(LF_vect2[dna(i)]['N']>LF_vect2[dna(i)][max]) max = 'N';
      if(LF_vect2[dna(i)]['T']>LF_vect2[dna(i)][max]) max = 'T';
      LF_max[i] = max;
    }

    for(uint64_t j = start; j <= end; ++j){

      if(bwt[j] != TERM){

        //if(std::find(high_freqs.begin(), high_freqs.end(), bwt[j]) == high_freqs.end()){	
        if(not high_freqs[ord(bwt[j])]){	

          if((int)(QUAL[j]-33) < quality_threshold){

            /*
             * If a low frequency base has also a low quality score we check if the symbol that precedes it is equal to the one that precedes one of the high frequency base.
             * In this case we replace the low frequency base with the high frequency base.
             */
            char freq1, freq2, new_base;
            int count = 0;

            //freq2 = std::max_element(LF_vect[bwt[j]].begin(), LF_vect[bwt[j]].end())-LF_vect[bwt[j]].begin();
            freq2 = LF_max[ord(bwt[j])];
            
            //for(auto it : high_freqs){
            for(int i=0; i<5; i++){
              if(high_freqs[i]){

                //freq1 = std::max_element(LF_vect[it].begin(), LF_vect[it].end())-LF_vect[it].begin();
                //freq1 = std::max_element(LF_vect[dna(i)].begin(), LF_vect[dna(i)].end())-LF_vect[dna(i)].begin();
                freq1 = LF_max[i]; 

                if(freq1 == freq2){
                  //new_base = it;
                  new_base = dna(i);
                  count++;
                } 
              }
            }

            if(count == 1){
              //FELIPE
              rankbv_setbit(rbv,j);
              //BWT_MOD[j] = new_base;
              BWT_MOD.push_back(new_base);
              modified++;
            }
          }
        }
        //FELIPE
        /*
        if(newqs != QUAL[j]) cout<<QUAL[j]<<" --> "<<newqs<<endl;
        else cout<<QUAL[j]<<endl;
        */
        QUAL[j] = newqs;
      }
    }
  }

  //reset temporary vector that stores frequencies in the cluster
  /*
  freqs['A'] = 0;
  freqs['C'] = 0;
  freqs['G'] = 0;
  freqs['T'] = 0;
  freqs['N'] = 0;
  high_freqs.clear();
  for(int i=0; i<5; i++){ //'A', 'C', 'G', 'T', 'N'
    freqs[i]=0;
    high_freqs[i]=false;
  }
  */    
  return;
}

/*
 * END PROCEDURE TO ANALYZE A CLUSTER
 */


/*
 * PROCEDURE run NAVIGATES suffix tree, and computes LCP minima, EXECUTES process_cluster for each detected cluster.
 */
void run(){
    
  ofstream out_file = ofstream(output);
  
  //read base qualities (permuted according to the BWT) into QUAL
  {
    ifstream f(input_qual); //taking file as inputstream
    if(f) {
      ostringstream ss;
      ss << f.rdbuf();
      QUAL = ss.str();
	    #if M==2
	      def_qs();
  	  #endif
    }
  }
  
  uint64_t n = bwt.size();

  //read BWT bases into the string BWT_MOD
  {
    /*
    ifstream f(input_dna); //taking file as inputstream
    if(f) {
        ostringstream ss;
        ss << f.rdbuf();
        BWT_MOD = ss.str();
    }
    */
    //The BWT is already in memory
    /*
    BWT_MOD = new char[n];
    for(uint64_t i=0; i<n; i++) 
      BWT_MOD[i]=bwt[i];
    */
    rbv = rankbv_create(n,2);
  }
  
  uint64_t begin = 0;//begin position
  
  uint64_t clust_len=0;
  bool cluster_open=false;
  
  int perc = -1;
  int last_perc = -1;
  
  uint64_t clust_size = 0; //cumulative cluster size
  
  #if DEBUG
    //used only to compute and visualize cluster statistics
    uint64_t MAX_CLUST_LEN = 200;
    //auto CLUST_SIZES = vector<uint64_t>(MAX_CLUST_LEN+1,0);
    uint64_t CLUST_SIZES[MAX_CLUST_LEN+1]{0};
  #endif
  
  //procedure that identifies clusters by looking at LCP_threshold and LCP_minima
  for(uint64_t i=0;i<n;++i){
      
    if(LCP_threshold[i] and not LCP_minima[i]){
          
      if(cluster_open){//extend current cluster
        clust_len++;
      }else{//open new cluster
        cluster_open=true;
        clust_len=1;
        begin=i;
      }
          
      }else{
          
        if(cluster_open){//close current cluster
          clust_size += clust_len;
          #if DEBUG
            if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]+=clust_len;
          #endif
          process_cluster(begin, i);//position i included
        }
        cluster_open=false;
        clust_len = 0;
      }
      
      perc = (100*i)/n;
      
      if(perc > last_perc){
        #if DEBUG
          if(verbose) cout << perc << "%. "<<endl;
        #endif
        last_perc = perc;
      }
  }
  
  cout << endl << "Done." << endl;

  //FELIPE
  rankbv_build(rbv);
  //cout << "rank = "<<rankbv_rank1(rbv,n)<<endl;
  
  //print clusters statistics (i.e. number of bases that fall inside each cluster of a fixed size)
  #if DEBUG
    uint64_t scale = *max_element(CLUST_SIZES.begin(), CLUST_SIZES.end());
    for(int i=0;i<=MAX_CLUST_LEN;++i){
      cout << i << ( i < 10 ? "   " : (i<100 ? "  " : " "));
      for(uint64_t j = 0; j < (100*CLUST_SIZES[i])/scale; ++j) cout << "-";
      cout << " " << CLUST_SIZES[i] << endl;
    }
  #endif
}

/*
 * END run
 */

/*
 * PROCEDURE invert INVERT BWT AND WRITE A NEW FASTQ FILE.
 *
 * If the input file consisted of reads and their reverse, we must define the behaviour.
 */
void invert(){

  /*
  ofstream out(output);
  ifstream in(original_fastq);
  */

  FILE *f_in = fopen(original_fastq.c_str(), "r");
  if(!f_in) perror("invert");

  FILE *f_out = fopen(output.c_str(), "w");
  if(!f_in) perror("invert");

  char header[]="@\n";
  char plus[]="+\n";
  size_t len = 0;

  char BASES[LONGEST]{0};
  char QS[LONGEST]{0};

  BASES[LONGEST-1] = '\n';
  QS[LONGEST-1] = '\n';


  //number of reads in the file
  //uint64_t N = bwt.rank(bwt.size(),TERM);
  uint64_t N = bwt.get_number_of_strings();
  
  string line;
  
  #if REVC == 0
    for(uint64_t i = 0;i < N;++i){//for each read 
         
      //string bases;
      //string qualities;
      
      int nbases = LONGEST-1;
      uint64_t j = i;//bwt[j] = current read character

      while(bwt[j] != TERM){
          
        //bases.push_back(BWT_MOD[j]);
        //BASES[--nbases] = BWT_MOD[j];
        BASES[--nbases] = (rankbv_access(rbv, j)==1)?(BWT_MOD[rankbv_rank1(rbv,j)-1]):(bwt[j]);

        #if B==1
          QUAL[j] = illumina_8_level_binning((int)QUAL[j]-33);
        #endif
        //qualities.push_back(QUAL[j]);
        QS[nbases] = QUAL[j];

        j = bwt.LF(j); //backward search
      }
      
      /*
      std::reverse(bases.begin(),bases.end());
      std::reverse(qualities.begin(),qualities.end());
      */

      #if DEBUG
        //for(auto q:qualities) statistics_qual_after[q-33]++;
        for(int i=LONGEST-nbases; i<LONGEST; i++) statistics_qual_after[q-33]++;
      #endif

      //write output FASTQ file
      char *buf = NULL;
    
      if(not ignore_headers){
        /*
        std::getline(in, line);//get read ID from the original FASTQ (headers)
        out << line << endl; //headers
        */
        ssize_t size = getline(&buf, &len, f_in); // @'s line
        fwrite(buf, sizeof(char), size, f_out);
      }
      else{
        /*
        out<<"@"<<endl; //default header
        */
        fwrite(header, sizeof(char), 2, f_out);
      }

      /*
      out << bases << endl; //bases
      out << "+" << endl;
      out << qualities << endl; //qs
      */

      //read input FASTQ file
      /*
      std::getline(in, line);//bases
      std::getline(in, line);//+
      std::getline(in, line);//qs
      */

      getline(&buf, &len, f_in); // bases 
      getline(&buf, &len, f_in); // +'s line
      getline(&buf, &len, f_in); // @'s line

      fwrite(&BASES[nbases], sizeof(char), LONGEST-nbases, f_out);
      fwrite(plus, sizeof(char), 2, f_out);
      fwrite(&QS[nbases], sizeof(char), LONGEST-nbases, f_out);

      free(buf);
      
      #if DEBUG
        for(auto q:line) statistics_qual_before[q-33]++;
      #endif
    }

  fclose(f_in);
  fclose(f_out);
  #else  
    //TODO
    /*
    for(uint64_t i = 0;i < (revc?N/2:N);++i){//for each read (if revc=true, only first half of reads)
         
      string bases;
      string qualities;
      
      uint64_t j = i;//bwt[j] = current read character
      
      while(bwt[j] != TERM){
          
        bases.push_back(BWT_MOD[j]);
        #if B==1
          QUAL[j] = illumina_8_level_binning((int)QUAL[j]-33);
        #endif
        qualities.push_back(QUAL[j]);
        j = bwt.LF(j); //backward search
      }
      
      std::reverse(bases.begin(),bases.end());
      std::reverse(qualities.begin(),qualities.end());
      
      //if second half of reads is the reverse complement of first half, combine the two
      if(revc){ 
        string bases_rc;
        string qualities_rc;
        
        j = i + N/2;//index of the corresponding read in the second half of the file
        
        while(bwt[j] != TERM){
          bases_rc.push_back(complement(BWT_MOD[j]));
          qualities_rc.push_back(QUAL[j]);
          j = bwt.LF(j);
        }
        
        if(bases_rc.length() != bases.length()){
          cout << "Error: second half of reads is not the reverse complement of first half. " << endl <<
          "found pair with different lengths (" << bases_rc.length() << "/" << bases.length() << ")" << endl;
          exit(0);
        }
        
        // DEFINE HOW TO COMBINE 
        for(int k=0;k<bases.length();++k){
          if (qualities[k] == qualities_rc[k]){
            if(bases[k] != bases_rc[k]){
              //??
            }
          }
          else{
            if(qualities[k] < qualities_rc[k]){
              bases[k] = bases_rc[k];
            }   
          }
        }//end-for  
      } //end if
      
      #if DEBUG
        for(auto q:qualities) statistics_qual_after[q-33]++;
      #endif

      std::getline(in, line);//get read ID from the original FASTQ (headers)
      
      //write output FASTQ file
      //cout << line << endl;

      out << line << endl; //headers
      out << bases << endl; //bases
      out << "+" << endl;
      out << qualities << endl; //qs
      
      //read input FASTQ file
      std::getline(in, line);//bases
      std::getline(in, line);//+
      std::getline(in, line);//qs
      
      #if DEBUG
        for(auto q:line) statistics_qual_before[q-33]++;
      #endif
    }//end for      
    */
  #endif
}
/*
 * END invert
 */

/*
 * PROCEDURES TO DEBUG
 */

// for each bwt position, the read coordinate, bwt base, modified base, and modified quality score
void print_info(){
    
  if(debug){
      
    //number of reads
    //uint64_t N = bwt.rank(bwt.size(),TERM);
    uint64_t N = bwt.get_number_of_strings();
    
    read_info = vector<string>(bwt.size());
    
    #if REVC
      for(uint64_t i = 0;i < (revc?N/2:N);++i){//for each read (if RC=true, only first half of reads)
    #else
      for(uint64_t i = 0;i < N;++i){//for each read (if RC=true, only first half of reads)
    #endif

      #if REVC
        for(int x=0;x<(revc?2:1);++x){
      #else
        int x=0;
        {
      #endif
        
        uint64_t j = i + x*(N/2);
        uint64_t off=0;//offset from the end of the read
        
        while(bwt[j] != TERM){
          read_info[j] = read_ids[i].substr(0,max_id_len);
          read_info[j].append(string("\t"));
          read_info[j].append(to_string(off));
          j = bwt.LF(j);
          off++;
        }
      }
    }
    
    cout << "ID\tposition\toriginal\tmodified\tmodified.quality\tLCP>=K\tminimum?" << endl;
    for(uint64_t i=0;i<bwt.size();++i){
      //cout << read_info[i] << "\t" << bwt[i] << "\t" << BWT_MOD[i] << "\t" << QUAL[i] << "\t" << (LCP_threshold[i]?"+\t":"\t") << (LCP_minima[i]?"*":"") << endl;
      cout << read_info[i] << "\t" << bwt[i] << "\t" << QUAL[i] << "\t" << (LCP_threshold[i]?"+\t":"\t") << (LCP_minima[i]?"*":"") << endl;
    }
      
  }
    
}
//load read IDs
void load_IDs(){
    
  ifstream in(original_fastq);

  string line;
  while (getline(in, line)){
    read_ids.push_back(line.substr(1));
    std::getline(in, line);//bases
    std::getline(in, line);//+
    std::getline(in, line);//qs
  }
}

/*
 * END PROCEDURES TO DEBUG
 */

int main(int argc, char** argv){
    
  srand(time(NULL));

  
  if(argc < 3) help();
  
  int opt;
  while ((opt = getopt(argc, argv, "he:q:o:f:k:m:t:rDvH")) != -1){
    switch (opt){
      case 'h':
        help();
        break;
      case 'e':
        input_dna = string(optarg);
        break;
      case 'q':
        input_qual = string(optarg);
        break;
      case 'o':
        output = string(optarg);
        break;
      case 'f':
        original_fastq = string(optarg);
        break;
      case 'k':
        K = atoi(optarg);
        break;
      case 'm':
        m = atoi(optarg);
        break;
      case 't':
        TERM = atoi(optarg);
        break;
      #if REVC
        case 'r':
          revc=true;
          break;
      #endif
      case 'D':
        debug=true;
        break;
      case 'H':
        ignore_headers=true;
        break;
      case 'v':
        verbose=true;
        break;
      default:
        help();
        return -1;
    }
  }

  K = K == 0 ? K_def : K;
  m = m == 0 ? m_def : m;

  if( input_dna.compare("")==0 or
      input_qual.compare("")==0 or
      output.compare("")==0 or
      original_fastq.compare("")==0
    ) help();

  if(not file_exists(input_dna)){
    cout << "Error: could not find file " << input_dna << "." << endl << endl;
    help();
  }

  if(not file_exists(input_qual)){
    cout << "Error: could not find file " << input_qual << endl << endl;
    help();
  }


  cout << "This is FASTQcompression." << endl;
  cout << "\tK: " << K << endl;
  cout << "Output fastq file: " << output << endl;
  cout << endl;

  cout << "Phase 1/4: loading and indexing eBWT ... " << flush;
  bwt = dna_bwt_n_t(input_dna,TERM);
  cout << "done." << endl;


  //number of reads in the file
  //uint64_t N = bwt.rank(bwt.size(),TERM);	
  uint64_t N = bwt.get_number_of_strings();
  cout << "Number of reads: " << N << endl;

  //detects clusters through local LCP minima
  detect_minima();

  //start procedure run			
  run();
  cout << "end run" << endl;
  //invert BWT
  invert();
  cout << "end invert" << endl;


  cout << clusters_size << " (" << (double(100*clusters_size)/bwt.size()) <<  "%) bases fall inside a cluster" << endl;
  cout << "done. " << modified << "/" << bwt.size() << " bases have been modified (" << 100*double(modified)/bwt.size() << "% of all bases and " <<
    100*double(modified)/clusters_size << "% of bases inside clusters)." << endl;

  if(debug)
    load_IDs();

  #if DEBUG
     cout << "Cumulative distribution of base qualities before: " << endl;
     uint64_t sum_tot = 0;
     for(auto x : statistics_qual_before) sum_tot+=x;

     uint64_t sum = 0;

     for(int i=0;i<50;++i){
     sum += statistics_qual_before[i];
     cout << i << "\t" << statistics_qual_before[i] << "\t" << double(sum)/sum_tot << endl;
     }
     cout << endl;

     cout << "Cumulative distribution of base qualities after: " << endl;
     sum_tot = 0;
     for(auto x : statistics_qual_after) sum_tot+=x;

     sum = 0;

     for(int i=0;i<50;++i){
     sum += statistics_qual_after[i];
     cout << i << "\t" << double(sum)/sum_tot << endl;
     }

     cout << endl;
  #endif

  if(debug)
    print_info();

  //delete[] LCP_minima;
  //delete[] LCP_threshold;
  //delete[] BWT_MOD;

  rankbv_free(rbv);

return 0;
}
