#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;


void readBL62(double* bl62[], char* matrix, double beta, short* letterToIndexMap);
double computeK1(double *bl62[]);

int main(int argc, char* argv[]) {
    //run with: exec.out bl62.txt 0.01
    double *bl62[20], S1;
    short indexMap[150]; //plug in letter, get index in bl62
    string beta(argv[2]);
    readBL62(bl62, argv[1], stod(beta), indexMap);
    S1 = computeK1(bl62);
} //main()

void readBL62(double* bl62[], char* matrix, double beta, short* indexMap) {
    //reads in BL62 from file, store in the values raised to power beta
    //updates indexMap to map letters to matrix locations
    //  can retrieve bl62[A][R] with:
    //      bl62[indexMap[A]][indexMap[R]]
    char letter;
    ifstream inf(matrix);

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
        } //read every value in the line, store the value^beta in matrix
    } //get every line, get rid of beginning letter

 } //readBL62() 

double computeK1(double* bl62[]) {
 
} //computeK1()

