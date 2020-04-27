// Strassen source: https://martin-thoma.com/strassen-algorithm-in-python-java-cpp/
// Sparse source: https://leetcode.com/accounts/login/?next=/problems/sparse-matrix-multiplication/

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cmath>

// Set LEAF_SIZE to 1 if you want to the pure strassen algorithm
// otherwise, the ikj-algorithm will be applied when the split
// matrices are as small as LEAF_SIZE x LEAF_SIZE
int leafsize = 16;
double total_time = 0;

using namespace std;

/*
 * Implementation of the strassen algorithm, similar to 
 * http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
 */

ostream& operator<<(ostream& os, const vector<vector<int>>& m)
{
    for (const auto& x : m) {
        for (const auto& y : x) {
            os << y << " ";
        }
        os << endl;
    }
    return os;
}


void strassen(vector< vector<int> > A, 
              vector< vector<int> > B, 
              vector< vector<int> > &C);

unsigned int nextPowerOfTwo(int n);

void strassenR(vector< vector<int> > &A, 
              vector< vector<int> > &B, 
              vector< vector<int> > &C, 
              int tam);

void sum(vector< vector<int> > &A, 
         vector< vector<int> > &B, 
         vector< vector<int> > &C, int tam);

void subtract(vector< vector<int> > &A, 
              vector< vector<int> > &B, 
              vector< vector<int> > &C, int tam);

void ikjalgorithm(vector< vector<int> > A, 
                  vector< vector<int> > B,
                  vector< vector<int> > &C, int n);

vector<vector<int>> read(string filename);

void sparseMultiply(vector< vector<int> > &A, 
                    vector< vector<int> > &B, 
                    vector< vector<int> > &C);

vector<vector<int>> transpose(const vector<vector<int>>& m) {
    int rows = m.size(), columns = m[0].size();
    vector<vector<int>> transposed_matrix(columns, vector<int>(rows, 0));

    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < rows; j++) {
            transposed_matrix[i][j] = m[j][i];
        }
    }

    return transposed_matrix;
}

int isSquare(const vector<vector<int>>& m) {
    if (m.size() == m[0].size()) {
        return 1;
    }
    else {
        return 0;
    }
}


int main (int argc, char* argv[]) {
    if (argc != 2) {
      cout << "Usage: ./a.out <file_name>\n";
      return -1;
    } 

    string filename = argv[1];
    vector<vector<int>> A = read(filename), B;
    
    if (!isSquare(A)) {
        B = transpose(A);
    }
    else {
        B = A;
    }

    vector<vector<int>> C(A.size(), vector<int>(B[0].size(), 0));
    vector<vector<int>> D = C, E = C;


    strassen(A, B, C);
    
    sparseMultiply(A, B, D);


    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < B[0].size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                E[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Regular: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0 << endl;
    

    
    return 0;
}


void ikjalgorithm(vector< vector<int> > A, 
                                   vector< vector<int> > B,
                                   vector< vector<int> > &C, int n) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void strassenR(vector< vector<int> >& A, 
              vector< vector<int> >& B, 
              vector< vector<int> > &C, int tam) {

    if (tam <= leafsize) {
        ikjalgorithm(A, B, C, tam);
        return;
    }

    // other cases are treated here:
    else {
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        int newTam = tam/2;
        vector<int> inner (newTam);
        vector< vector<int> > 
            a11(newTam,inner), a12(newTam,inner), a21(newTam,inner), a22(newTam,inner),
            b11(newTam,inner), b12(newTam,inner), b21(newTam,inner), b22(newTam,inner),
              c11(newTam,inner), c12(newTam,inner), c21(newTam,inner), c22(newTam,inner),
            p1(newTam,inner), p2(newTam,inner), p3(newTam,inner), p4(newTam,inner), 
            p5(newTam,inner), p6(newTam,inner), p7(newTam,inner),
            aResult(newTam,inner), bResult(newTam,inner);

        int i, j;

        //dividing the matrices in 4 sub-matrices:
        for (i = 0; i < newTam; i++) {
            for (j = 0; j < newTam; j++) {
                a11[i][j] = A[i][j];
                a12[i][j] = A[i][j + newTam];
                a21[i][j] = A[i + newTam][j];
                a22[i][j] = A[i + newTam][j + newTam];

                b11[i][j] = B[i][j];
                b12[i][j] = B[i][j + newTam];
                b21[i][j] = B[i + newTam][j];
                b22[i][j] = B[i + newTam][j + newTam];
            }
        }

        // Calculating p1 to p7:
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        total_time += (chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0);

        sum(a11, a22, aResult, newTam); // a11 + a22
        sum(b11, b22, bResult, newTam); // b11 + b22
        strassenR(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)

        sum(a21, a22, aResult, newTam); // a21 + a22
        strassenR(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)

        subtract(b12, b22, bResult, newTam); // b12 - b22
        strassenR(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)

        subtract(b21, b11, bResult, newTam); // b21 - b11
        strassenR(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)

        sum(a11, a12, aResult, newTam); // a11 + a12
        strassenR(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)   

        subtract(a21, a11, aResult, newTam); // a21 - a11
        sum(b11, b12, bResult, newTam); // b11 + b12
        strassenR(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)

        subtract(a12, a22, aResult, newTam); // a12 - a22
        sum(b21, b22, bResult, newTam); // b21 + b22
        strassenR(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)

        // calculating c21, c21, c11 e c22:

        sum(p3, p5, c12, newTam); // c12 = p3 + p5
        sum(p2, p4, c21, newTam); // c21 = p2 + p4

        sum(p1, p4, aResult, newTam); // p1 + p4
        sum(aResult, p7, bResult, newTam); // p1 + p4 + p7
        subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7

        sum(p1, p3, aResult, newTam); // p1 + p3
        sum(aResult, p6, bResult, newTam); // p1 + p3 + p6
        subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6

        begin = chrono::steady_clock::now();
        // Grouping the results obtained in a single matrix:
        for (i = 0; i < newTam ; i++) {
            for (j = 0 ; j < newTam ; j++) {
                C[i][j] = c11[i][j];
                C[i][j + newTam] = c12[i][j];
                C[i + newTam][j] = c21[i][j];
                C[i + newTam][j + newTam] = c22[i][j];
            }
        }
        end = chrono::steady_clock::now();
        total_time += (chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0);
    }
}


unsigned int nextPowerOfTwo(int n) {
    return pow(2, int(ceil(log2(n))));
}


/* void strassen(vector< vector<int> > &A, 
              vector< vector<int> > &B, 
              vector< vector<int> > &C, unsigned int n) {
    //unsigned int n = tam;
    unsigned int m = nextPowerOfTwo(n);
    vector<int> inner(m);
    vector< vector<int> > APrep(m, inner), BPrep(m, inner), CPrep(m, inner);

    for(unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) {
            APrep[i][j] = A[i][j];
            BPrep[i][j] = B[i][j];
        }
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    strassenR(APrep, BPrep, CPrep, m);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Strassen: " << (chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0) - total_time << endl;
    

    for(unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) {
            C[i][j] = CPrep[i][j];
        }
    }
} */

void strassen(vector<vector<int>> mtx, vector<vector<int>> mtx2, vector<vector<int>>& newMtx) {
    int rows = mtx.size();
    int cols = mtx[0].size();
    int rows2 = mtx2.size();
    int cols2 = mtx2[0].size();
    double temp;

    // TODO:: Make it so it is a square matrix with 0's added
    if (rows < cols) {
        rows = cols;
    }
    if (rows2 < cols2) {
        rows2 = cols2;
    }
    if (rows < rows2) {
        rows = rows2;
    }

    // Calculate the next highest power of 2
    temp = ceil(log2(rows));
    rows = pow(2, temp);

    // Resizing and padding the new matrices for Strassen's
    mtx.resize(rows);
    mtx2.resize(rows);
    newMtx.resize(rows);
    for(int i = 0; i < rows; ++i) {
        mtx[i].resize(rows);
        mtx2[i].resize(rows);
        newMtx[i].resize(rows);
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    strassenR(mtx, mtx2, newMtx, rows);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Strassen: " << (chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0) - total_time << endl;

}

void sum(vector< vector<int> > &A, 
         vector< vector<int> > &B, 
         vector< vector<int> > &C, int tam) {
    int i, j;

    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}


void subtract(vector< vector<int> > &A, 
              vector< vector<int> > &B, 
              vector< vector<int> > &C, int tam) {
    int i, j;

    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }   
}

int getMatrixSize(string filename) {
    string line;
    ifstream infile;
    infile.open (filename.c_str());
    getline(infile, line);
    return count(line.begin(), line.end(), '\t') + 1;
}


vector<vector<int>> read(string filename) {
    ifstream fin;
    fin.open(filename.c_str());

    int row, col;
    
    fin >> row; fin >> col;

    vector<vector<int>> m(row, vector<int>(col, 0));

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            fin >> m[i][j];
        }
    }

    fin.close();

    return m;

}


void sparseMultiply(vector< vector<int> > &A, vector< vector<int> > &B, vector< vector<int> > &C) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

     for (int i = 0; i < C.size(); i++){
        for (int k = 0; k < A[0].size(); k++){
            if (A[i][k] != 0) {
                for (int j =0; j < C[0].size(); j++){
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    cout << "Sparse: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0 << endl;
}
