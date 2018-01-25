#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;


void readBL62(double* bl62[], string matrix, double beta, short* letterToIndexMap);
double computeK1(double *bl62[]);
void readInput(string matrix, string b, string seq1, string seq2);
string readSequence(string s);
double computeK1(double* bl62[], short* indexMap, string seq1, string seq2);

int main(int argc, char* argv[]) {
    //run with: exec.out bl62.txt Beta=0.01 seq1 seq2
	readInput(argv[1], argv[2], argv[3], argv[4]);// reads in arguments passed to program
   	return 0;
} //main()

void readInput(string matrix, string b, string seq1, string seq2){
	double *bl62[20];
	short indexMap[150]; //plug in letter, get index in bl62
	string beta(b);
	double sum1; //sum of first kernel

	readBL62(bl62, matrix, stod(beta), indexMap);//stod error
	seq1 = readSequence(seq1);
	cout<< seq1<< endl;
	seq2 = readSequence(seq2);
	cout<< seq2<< endl;

	sum1 = computeK1(bl62,indexMap, seq1, seq2);//computes the first kernel
	cout<< sum1;
}//readInput()

string readSequence(string s){
	//reads sequence from .fa file and returns string
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
    //  can retrieve bl62[A][R] with:
    //      bl62[indexMap[A]][indexMap[R]]
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
            cout << bl62[i][j] << ","; //debug
        } //read every value in the line, store the value^beta in matrix
        cout << "EOL" << endl; //debug (End of Line)
    } //get every line, get rid of beginning letter

 } //readBL62() 

double computeK1(double* bl62[], short* indexMap, string seq1, string seq2){
 //computes the first Kernel by comparing a single letter from
 //seq1 and seq2 and finding the bl62 value of the amino acids
 double sum;// sum of k1 values for all the single letter pairs
 unsigned int s1_size = seq1.length();
 unsigned int s2_size = seq2.length(); 
 for(unsigned int i =0; i < s1_size; i++){//warning message
	 for(unsigned int j = 0; j < s2_size; j++){
		sum= sum + bl62[indexMap[seq1[i]]][indexMap[seq2[j]]];
	 }
 }

	return sum;
} //computeK1()

