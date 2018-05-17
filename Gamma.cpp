#include "MatrixMult.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <limits>
struct Point
{
        double X, Y, Dose;
};


Matrix ReadMatrix(std::string Path) noexcept
{
        std::ifstream Input{Path};
        size_t Dim[2];
        double Coor[2];
        double VoxSize[2];
        for (size_t i=0;i<2;++i) Input >> Dim[i] >> Coor[i] >> VoxSize[i];
        Matrix Image(Dim[0],Dim[1]);
        Image.SetXCorner(Coor[0]);
        Image.SetYCorner(Coor[1]);
        Image.SetSizeX(VoxSize[0]);
        Image.SetSizeY(VoxSize[1]);
        for (size_t i=0;i<Dim[0];++i)
        {
                for (size_t j=0;j<Dim[1];++j) Input >> Image(i,j);
        }
        return Image;
}

double GetMax(const Matrix& M) noexcept
{
        size_t Rows=M.GetRows();
        size_t Columns=M.GetColumns();
        double Max=0;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0;j<Columns;++j)
                {
                        if (M(i,j)>Max) Max=M(i,j);
                }
        }
        return Max;
}

double GetMin(const Matrix& M) noexcept
{
        size_t Rows=M.GetRows();
        size_t Columns=M.GetColumns();
        double Min=DBL_MAX;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0;j<Columns;++j)
                {
                        if (std::isnan(M(i,j))) continue;
                        if (M(i,j)<Min) Min=M(i,j);
                }
        }
        return Min;
}

std::vector<size_t> SortRows(const std::vector<Point>& E, const double* const R) noexcept
{
        using std::vector;
        using std::sort;
        using std::pow;
        vector<double> Dist;
        size_t Size=E.size();
        for (size_t i=0;i<Size;++i) Dist.push_back(pow(R[0]-E[i].X,2)+pow(R[1]-E[i].Y,2)+pow(R[2]-E[i].Dose,2));
        vector<double> TmpDist=Dist;
        vector<size_t> ReturnVect;
        sort(Dist.begin(),Dist.end());
        for (size_t i=0;i<Size;++i) 
        {
                for (size_t j=0;j<Size;++j)
                {
                        if (Dist[i]==TmpDist[j])
                        {
                                ReturnVect.push_back(j);
                                break;
                        }
                }
        }
        return ReturnVect;
}

std::vector<Point> CreateSquare(const Matrix& Eva, const size_t Row, const size_t Column, const double D, const double d) noexcept
{
        std::vector<Point> Square;
        Point EntryPoint;

        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Dose=Eva(Row,Column)/D;
        Square.push_back(EntryPoint);

        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Dose=Eva(Row,Column+1)/D;
        Square.push_back(EntryPoint);

        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Dose=Eva(Row+1,Column+1)/D;
        Square.push_back(EntryPoint);
        
        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Dose=Eva(Row+1,Column)/D;
        Square.push_back(EntryPoint);

        return Square;
}

Matrix SetUpP(std::vector<Point> SimplexVertices, const double* const R) noexcept
{
        size_t Size=SimplexVertices.size();//SimplexVertices.size();
        Matrix P(3,1);
        P(0,0)=R[0]-SimplexVertices[Size-1].X;
        P(1,0)=R[1]-SimplexVertices[Size-1].Y;
        P(2,0)=R[2]-SimplexVertices[Size-1].Dose;

        return P;
}

Matrix SetUpV(std::vector<Point> SimplexVertices) noexcept
{
        size_t Size=SimplexVertices.size();
        Matrix V(3,Size-1);
        for (size_t VCol=0;VCol<(Size-1);++VCol)
        {
                V(0,VCol)=SimplexVertices[VCol].X-SimplexVertices[Size-1].X;
                V(1,VCol)=SimplexVertices[VCol].Y-SimplexVertices[Size-1].Y;
                V(2,VCol)=SimplexVertices[VCol].Dose-SimplexVertices[Size-1].Dose;
        }
        return V;
}


Matrix MLDivide(Matrix& A,Matrix& B) noexcept
{
       A.GaussianElimination(B);
       Matrix x = A.BackSubstitute(B);
       return x;
}

Matrix GetW(Matrix& VTV, Matrix& VT, const Matrix& P) noexcept
{
        Matrix TmpW=(MLDivide(VTV,VT))*P;
        Matrix W(TmpW.GetRows()+1,1);
        double Sum=0;
        for (size_t i=0;i<TmpW.GetRows();++i) 
        {
                W(i,0)=TmpW(i,0);
                Sum+=TmpW(i,0);
        }
        W(TmpW.GetRows(),0)=1-Sum;
        return W;
}

Point ComputeNearestPoint(const Matrix& W, const std::vector<Point>& SimplexVertices) noexcept
{
        Point NearestPoint;
        for (size_t i=0;i<W.GetRows();++i) 
        {
                NearestPoint.X+=W(i,0)*SimplexVertices[i].X;
                NearestPoint.Y+=W(i,0)*SimplexVertices[i].Y;
                NearestPoint.Dose+=W(i,0)*SimplexVertices[i].Dose;
        }
        return NearestPoint;
}

double GetGammaFromLine(const std::vector<Point>& Line, const Matrix& W, const double * const R) noexcept
{
        using std::pow;
        bool PositiveW=true;
        for (size_t i=0;i<W.GetRows();++i)
        {
                if (W(i,0)<0)
                {
                        PositiveW=false;
                        break;
                }
        }
        if (PositiveW) 
        {
                Point NP=ComputeNearestPoint(W,Line);
                return pow(R[0]-NP.X,2)+pow(R[1]-NP.Y,2)+pow(R[2]-NP.Dose,2);
        }
        else 
        {
                Point P1, P2;
                P1.X=Line[0].X;
                P1.Y=Line[0].Y;
                P1.Dose=Line[0].Dose;
                P2.X=Line[1].X;
                P2.Y=Line[1].Y;
                P2.Dose=Line[1].Dose;
                double GAMMA1=pow(R[0]-P1.X,2)+pow(R[1]-P1.Y,2)+pow(R[2]-P1.Dose,2);
                double GAMMA2=pow(R[0]-P2.X,2)+pow(R[1]-P2.Y,2)+pow(R[2]-P2.Dose,2);
                return GAMMA1<GAMMA2 ? GAMMA1 : GAMMA2;
        }

}

double CheckLine(const std::vector<Point> Line, const double* const R) noexcept
{
  
        Matrix P=SetUpP(Line,R);
        Matrix V=SetUpV(Line);
        Matrix VT=V.Transpose();
        Matrix VTV=VT*V;
        Matrix W=GetW(VTV,VT,P);
        return GetGammaFromLine(Line,W,R); 
}

double GetGammaFromTriangle(const std::vector<Point>& Triangle,const Matrix& W, const double * const R) noexcept
{

        bool PositiveW=true;
        for (size_t i=0;i<W.GetRows();++i)
        {
                if (W(i,0)<0)
                {
                        PositiveW=false;
                        break;
                }
        }
        double GAMMA;
        if (PositiveW) 
        {
                Point NP=ComputeNearestPoint(W,Triangle);
                GAMMA=pow(R[0]-NP.X,2)+pow(R[1]-NP.Y,2)+pow(R[2]-NP.Dose,2);
        }
        else 
        {
                std::vector<Point> NonNegW;
                Point Entry;
                std::vector<double> Gamma;
                for (size_t i=0;i<W.GetRows();++i)
                {
                        if (W(i,0)<0)
                        {
                                for (size_t j=0;j<W.GetRows();++j)
                                {
                                        if (i!=j)
                                        {
                                                Entry.X=Triangle[j].X;
                                                Entry.Y=Triangle[j].Y;
                                                Entry.Dose=Triangle[j].Dose;
                                                NonNegW.push_back(Entry);
                                        } 
                                }
                                GAMMA=CheckLine(NonNegW, R);
                                NonNegW.clear();
                        }
                } 
        }
        return GAMMA;
}

double CheckTriangles(const Matrix& Eva, const double D, const double d, const double* const R, const size_t Row, const size_t Column) noexcept
{
        using std::vector;
        using std::pow;
        
        vector<Point> Square=CreateSquare(Eva,Row,Column,D,d);
        vector<size_t> SortedRows=SortRows(Square,R);

        vector<Point> Triangle1, Triangle2;

        Triangle1.push_back(Square[SortedRows[0]]);
        Triangle1.push_back(Square[SortedRows[1]]);
        Triangle1.push_back(Square[SortedRows[2]]);

        Matrix P=SetUpP(Triangle1, R);
        Matrix V=SetUpV(Triangle1);
        Matrix VT=V.Transpose();
        Matrix VTV=VT*V;
        Matrix W=GetW(VTV,VT,P);
        double GAMMA1=GetGammaFromTriangle(Triangle1,W,R);
        
        
        Triangle2.push_back(Square[SortedRows[1]]);
        Triangle2.push_back(Square[SortedRows[3]]);
        Triangle2.push_back(Square[SortedRows[2]]);

        P=SetUpP(Triangle2, R);
        V=SetUpV(Triangle2);
        VT=V.Transpose();
        VTV=VT*V;
        W=GetW(VTV,VT,P);
        double GAMMA2=GetGammaFromTriangle(Triangle2,W,R);
        return GAMMA1<GAMMA2 ? GAMMA1 : GAMMA2;
}

void SetValue(Matrix& M, double Value, size_t Row, int Column) noexcept
{
        if (Column==-1)
        {
                for (size_t Col=0; Col<M.GetColumns();++Col) M(Row,Col)=Value;
        }
        else M(Row,Column)=Value;
}

Matrix Gamma2D(const Matrix& Ref, const Matrix& Eva, double Dose, const double d, const double DoseLim, const double SearchLim) noexcept
{
        using std::abs;
        using std::sqrt;
        Matrix Gamma(Ref.GetRows(),Ref.GetColumns(),std::numeric_limits<double>::quiet_NaN());
        Matrix GAMMA(Eva.GetRows(),Eva.GetColumns(),std::numeric_limits<double>::quiet_NaN());

        const double MaxDose=GetMax(Ref);
        const double D=((100+Dose)/100)*MaxDose-MaxDose;
        const double DL=((100+DoseLim)/100)*MaxDose-MaxDose;
        const double SL=SearchLim/d;

        double* Rr = new double[3];
        double Er;

        size_t RefRows=Ref.GetRows();
        size_t RefColumns=Ref.GetColumns();
        size_t EvaRows=Eva.GetRows();
        size_t EvaColumns=Eva.GetColumns();
        for (size_t RefRow=0; RefRow<RefRows;++RefRow)
        {
                Rr[1]=Ref.GetYPosition(RefRow)/d;
                for (size_t RefCol=0; RefCol<RefColumns;++RefCol)
                {
                        if (Ref(RefRow,RefCol)<DL) continue;
                        Rr[0]=Ref.GetXPosition(RefCol)/d;
                        Rr[2]=Ref(RefRow,RefCol)/D;
                        for (size_t EvaRow=0; EvaRow<(EvaRows-1);++EvaRow) 
                        {
                                Er=Eva.GetYPosition(EvaRow)/d;
                                if (abs(Er-Rr[1])>SL) 
                                {
                                        SetValue(GAMMA,std::numeric_limits<double>::quiet_NaN(),EvaRow,-1);
                                        continue;
                                }
                                for (size_t EvaCol=0;EvaCol<(EvaColumns-1);++EvaCol)
                                {
                                        Er=Eva.GetXPosition(EvaCol)/d;
                                        if (abs(Er-Rr[0])>SL) 
                                        {
                                                
                                                SetValue(GAMMA,std::numeric_limits<double>::quiet_NaN(),EvaRow,EvaCol);
                                                continue;
                                        }
                                        GAMMA(EvaRow,EvaCol)=CheckTriangles(Eva,D,d,Rr,EvaRow,EvaCol);                                         
                                }
                        }
                        Gamma(RefRow,RefCol)=sqrt(GetMin(GAMMA));
                }
        }
        return Gamma;
}
