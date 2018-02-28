/*
 * kernel.cpp Sequence Kernel Program 
 * By: Zaid Al Rakabi and Munir Nur
*/
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;


void readBL62(double* bl62[], string matrix, double beta, short* letterToIndexMap);
double computeKhat(double *bl62[], short* indexMap, string seq1, string seq2);
void readInput(string matrix, char* b, string seq1, string seq2);
string readSequence(string s);
double computeK(double* bl62[], short* indexMap, string seq1, string seq2);
double sumK(double* bl62[], short* indexMap, string seq1, string seq2, int k);


int main(int argc, char* argv[]) {
    //run with: exec.out bl62.txt Beta=0.01 seq1 seq2
	readInput(argv[1], argv[2], argv[3], argv[4]);// reads in arguments passed to program
   	return 0;
} //main()

void readInput(string matrix, char* b, string seq1, string seq2){
	double *bl62[20];
	short indexMap[150]; //plug in letter, get index in bl62
	char* beta(b);
	double khat; //normalized score for 2 protein sequences similarity
	double dist; //dist = 0 means the sequences are similar, dist = 1 means they are not similar

	readBL62(bl62, matrix, atof(beta), indexMap);
	seq1 = readSequence(seq1);
	seq2 = readSequence(seq2);
	
	clock_t begin= clock();

	khat = computeKhat(bl62,indexMap, seq1, seq2);//computes khat

	dist = sqrt(2*(1 - khat)); //computes distance of khat for score of similarity

	clock_t end= clock();
	cout<<"distance: " <<  dist << endl;//print distance
	cout << "CPU time: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;//print time
}//readInput()

string readSequence(string s){
	//reads sequence from file and returns string
	// of all the letters of the sequence. 
	string seqstr= "";
	string tempstr= "";
	ifstream fstream(s.c_str());
	fstream.ignore(500, '\n');//ignore first line of file
	while(getline(fstream, tempstr)){
		seqstr = seqstr + tempstr;//concatenate lines if file has multiple lines
	}

	return seqstr;
}//readSequence()

void readBL62(double* bl62[], string matrix, double beta, short* indexMap) {
    //reads in BL62 from file, store in the values raised to power beta
    //updates indexMap to map letters to matrix locations
    //can retrieve bl62[A][R] with:
    //bl62[indexMap[A]][indexMap[R]]
    unsigned char letter;
    ifstream inf(matrix.c_str());

    for (int i = 0; i < 20; i++)
        bl62[i] = new double[20]; //allocate matrix

    for (int i = 0; i < 20; i++) {
        inf >> letter;
        indexMap[letter] = i;
    } //turn letters into indices

    for (int i = 0; i < 20; i++) {
        inf >> letter;
        for (int j = 0; j < 20; j++) {
            inf >> bl62[i][j];
            bl62[i][j] = pow(bl62[i][j], beta);
         //test print:   cout << bl62[i][j] << ","; //debug
        } //read every value in the line, store the value^beta in matrix
       //test print: cout << "EOL" << endl; //debug (End of Line)
    } //get every line, get rid of beginning letter

 } //readBL62()

double computeKhat(double* bl62[], short* indexMap, string seq1, string seq2){
	//computes a normalized score for the K value using the formula for Khat
	double kST, kSS, kTT, kH, d;

	kST= computeK(bl62, indexMap, seq1, seq2);
	kSS= computeK(bl62, indexMap, seq1, seq1);
	kTT= computeK(bl62, indexMap, seq2, seq2);
	kSS= kSS * kTT;

	kH= (kST/sqrt(kSS));//formula: k(s1,s2)/sqrt(k(s1,s1)*k(s2,s2))
	
	return kH;
}//computeKhat

double computeK(double* bl62[], short* indexMap, string seq1, string seq2){
 //computes all the Kernels by comparing a single letter from
 //seq1 and seq2 and finding the bl62 value of the amino acids
 double sum;// sum of k values for all the possible pairs
 int s1_size = seq1.length();
 int s2_size = seq2.length();
 int s_max = min(s1_size, s2_size);

 double total_K = 0;
 for(int k = 1; k<= s_max; k++){//find sum for all k values up to s_max
	 sum = sumK(bl62, indexMap, seq1, seq2, k);
	 total_K += sum; 
 }
  return total_K;
}//computeK()

double sumK(double* bl62[], short* indexMap, string seq1, string seq2, int k){
 //computes the sum of kth Kernel by comparing a k elements of the sequence
 //using seq1 and seq2 and finding the bl62 value of the amino acids
 double sum;// sum of amino acid values for all the k letter pairs
 unsigned int s1_size = seq1.length();
 unsigned int s2_size = seq2.length();

 for(unsigned int i = 0; i < (s1_size - k); i++){//warning message

	 for(unsigned int j = 0; j < (s2_size -k); j++){
         double temp = 1;

         for (int l = 0; l < k; l++) {
             temp = temp * bl62[indexMap[seq1[i+ l]]][indexMap[seq2[j+ l]]];
         }	
         sum = sum + temp;
     }

 }

	return sum;
} //sumK()

