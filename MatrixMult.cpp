#include "MatrixMult.h"
#include <iomanip>
Matrix::Matrix(size_t NoRows, size_t NoColumns, double Val)
        :Rows{NoRows}, Columns{NoColumns}, Mat{new double[NoRows*NoColumns]}
{
        for (size_t i=0;i<Rows*Columns;++i) Mat[i]=Val;
}


Matrix::Matrix(Matrix&& Rhs)
{
        Rows=Rhs.Rows;
        Columns=Rhs.Columns;
        Mat=Rhs.Mat;
        Rhs.Mat=nullptr;
        Rhs.Rows=0;
        Rhs.Columns=0;
}

Matrix& Matrix::operator=(Matrix&& Rhs)
{
        if (this !=&Rhs)
        {
                delete[] Mat;
                Mat = Rhs.Mat;
                Rows=Rhs.Rows;
                Columns=Rhs.Columns;

                Rhs.Mat=nullptr;
                Rhs.Rows=0;
                Rhs.Columns=0;
        }
        return *this;
}

void Matrix::PrintMatrix() const noexcept
{
        using std::cout;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0; j<Columns; ++j) cout << "  " << (*this)(i,j);
                cout << '\n';
        }
        cout << '\n';
}

Matrix Matrix::operator*(const Matrix& Rhs)
{
        using std::cout;
        if (Rows!=Rhs.GetColumns() && Columns!=Rhs.GetRows())
        {
                cout << "Dimensions of the matrices are incorrect!\n";
                exit (1);
        }
        Matrix ReturnMatrix(Rows,Rhs.GetColumns());
        for (size_t i=0;i<ReturnMatrix.GetRows();++i)
        {
                for (size_t j=0;j<ReturnMatrix.GetColumns();++j) 
                {
                        for (size_t k=0;k<Columns;++k) ReturnMatrix(i,j)+=(*this)(i,k)*Rhs(k,j);
                }
        }
        return ReturnMatrix;
}

Matrix Matrix::Transpose() const noexcept
{
        Matrix TempMatrix(Columns,Rows);
        for (size_t i=0;i<Columns;i++)
        {
                for (size_t j=0;j<Rows;++j) TempMatrix(i,j)=(*this)(j,i);
        }
        return TempMatrix;
}

void Matrix::SwapRows(size_t Row1, size_t Row2) noexcept
{
       if (Row1==Row2) return;
       double TempElement;
       for (size_t i=0;i<Columns;++i)
       {
               TempElement=(*this)(Row1,i); 
               (*this)(Row1,i)=(*this)(Row2,i);
               (*this)(Row2,i)=TempElement;
       }       
}

void Matrix::GaussianElimination(Matrix& B) noexcept
{
        double MaxEl;
        size_t MaxRow;
        size_t BColumns=B.GetColumns();
        for (size_t i=0;i<Rows;++i)
        {
                MaxEl=(*this)(i,i);
                MaxRow=i;
                for (size_t j=i+1;j<Rows;++j)
                {
                        if ((*this)(j,i)>MaxEl)
                        {
                                MaxEl=(*this)(j,i);
                                MaxRow=j;
                        }
                }
                SwapRows(i,MaxRow);
                B.SwapRows(i,MaxRow);
                double FirstElement;
                for (size_t j=i+1;j<Rows;++j)
                {
                        for (size_t k=0;k<BColumns;++k) B(j,k)=B(j,k)-B(i,k)*(*this)(j,i)/MaxEl;
                        FirstElement=(*this)(j,i);
                        for (size_t k=0;k<Columns;++k)
                        {
        
                                (*this)(j,k)=(*this)(j,k)-(*this)(i,k)*FirstElement/MaxEl;
                        }
                }       
        }
}

Matrix Matrix::BackSubstitute(Matrix& B) noexcept
{
        Matrix x(B.GetRows(),B.GetColumns());

        for (int i=(B.GetRows()-1);i>=0;--i) 
        {
                for (size_t j=0;j<B.GetColumns();++j)
                {
                        x(i,j)=B(i,j)/(*this)(i,i);
                        for (int k=i-1;k>=0;--k)
                        {
                                B(k,j)=B(k,j)-x(i,j)*(*this)(k,i);
                        }
                }
        }
        return x;
}

/*Matrix Matrix::MLDivide(Matrix& B) noexcept
{
       (*this).GaussianElimination(B);
       Matrix x = (*this).BackSubstitute(B);
       return x;
}*/


