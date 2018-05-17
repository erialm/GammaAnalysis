#include "MatrixMult.h"
#include <iostream>
//#include <chrono>
int main()
{
        using std::string;

        Matrix ReadMatrix(string); //forward declare
        Matrix Gamma2D(const Matrix&, const Matrix&, const double, const double, const double, const double);
        Matrix Ref=ReadMatrix("../32x32RefMatrix.txt"); //32x32
        Matrix Eva=ReadMatrix("../241x241EvaMatrix.txt"); //241x241
        double D=3;
        double d=2;
        double DoseLim=10;
        double SearchLim=4.5;
        //auto start = std::chrono::high_resolution_clock::now();
        Matrix Gamma=Gamma2D(Ref,Eva,D,d,DoseLim,SearchLim);
        //auto finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;
        //std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        Gamma.PrintMatrix();
        return 0;
}
