#include <iostream>
#include <fstream>
#include <string>
class Matrix
{
        private:
                double* Mat;
                size_t Rows, Columns;
                double XCorner, YCorner;
                double SizeX, SizeY;

        public:
                Matrix(size_t Rows, size_t Columns, double Val=0); //Constructor
                Matrix(const Matrix&) = delete; //No copy constructor
                Matrix (Matrix&&);    //Move constructor

                ~Matrix() {if (Mat!=nullptr) delete[] Mat;}      //Destructor
                
                size_t GetRows() const noexcept {return Rows;} 
                size_t GetColumns() const noexcept {return Columns;}
                double GetXPosition(size_t Column) const noexcept {return XCorner+Column*SizeX;}
                double GetYPosition(size_t Row) const noexcept {return YCorner+Row*SizeY;}

                void SetXCorner(double X) {XCorner = X;}
                void SetYCorner(double Y) {YCorner = Y;}
                void SetSizeX(double X) {SizeX = X;}
                void SetSizeY(double Y) {SizeY = Y;}
                 

                double  operator() (size_t Row, size_t Column) const {return Mat[Columns*Row + Column];}  
                double  operator() (size_t Index) const {return Mat[Index];}  
                double& operator() (size_t Row, size_t Column) {return Mat[Columns*Row + Column];}       

                Matrix& operator=(const Matrix&) = delete; //No copy assignment
                Matrix operator*(const Matrix&); //Matrix multiplication
                Matrix& operator=(Matrix&&);

                void PrintMatrix() const noexcept;
                Matrix Transpose() const noexcept;
                void SwapRows(size_t,size_t) noexcept;
                void GaussianElimination(Matrix&) noexcept;
                Matrix BackSubstitute(Matrix&) noexcept;
                //Matrix MLDivide(Matrix&) noexcept;
};



