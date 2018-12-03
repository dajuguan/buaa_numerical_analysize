//本文件是数值分析第二个作业定义的通用矩阵操作头文件matrix.h
//Copyright (c) 2018-2020 chenjunfeng <wwwcjf163@163.com>
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.

#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>
#include<iomanip>
#include<assert.h>

using namespace std;
#define TOL 1e-12
#define MAX_STEPS 300

template<typename T>
class Matrix;

//实例化double类型的矩阵
typedef Matrix<double> Mat;

template<typename T>
ostream &operator<< (ostream&out, const Matrix<T>&obj){  // 输出到屏
    for(int i=0;i<obj.Row();i++){
        out<<'\t';
        for(int j=0; j<obj.Col();j++){
            out<<setiosflags(ios::showpos|ios::fixed)<<setprecision(8)<<obj[i][j]<<"\t";
        }
        out<<endl;
    }
    return out;
}  

//矩阵相乘
template<typename T>
Matrix<T> operator*(const Matrix<T>&lm, const Matrix<T>&rm){
    assert(lm.col == rm.row);
    Matrix<T> ret(lm.row,rm.col);
    for(int i=0;i<lm.row;i++){
        for(int j=0;j<rm.col;j++){
            for(int k=0;k<lm.col;k++){
                ret.pmm[i*rm.col + j] += lm[i][k]*rm[k][j];
            }
        }
    }
    return ret;
}    

//矩阵左乘以常数字
template<typename T>
Matrix<T> operator*(const Matrix<T>&lm, T val){
    Matrix<T> ret(lm.row,lm.col);
    for(int i=0;i<ret.size;i++){
        ret.pmm[i] = val * lm.pmm[i];
    }
    return ret;
}

//除以常数
template<typename T>
Matrix<T> operator/(const Matrix<T>&lm, T val){
    T inverseval = 1/val;
    return lm*inverseval;
}

//矩阵相减
template<typename T>
Matrix<T> operator- (const Matrix<T>&lm,  const Matrix<T>&rm){
    assert(lm.row==rm.row && lm.col == rm.col);
    Matrix<T> ret(lm.row,lm.col);
    for(int i=0;i<lm.size;i++){
        ret.pmm[i] = lm.pmm[i] - rm.pmm[i];
    }
    return ret;
}

//矩阵相加
template<typename T>
Matrix<T> operator+ (const Matrix<T>&lm,  const Matrix<T>&rm){
    assert(lm.row==rm.row && lm.col == rm.col);
    Matrix<T> ret(lm.row,lm.col);
    for(int i=0;i<lm.size;i++){
        ret.pmm[i] = lm.pmm[i] + rm.pmm[i];
    }
    return ret;
}

//复数类
class Complex{
    public:
        double real,imag;
        Complex():real(0),imag(0){}
        Complex(double r,double i):real(r),imag(i){}
        Complex& operator=(const Complex&rhs){
            this->real = rhs.real;
            this->imag = rhs.imag;
            return *this;
        }
        Complex& operator=(const double&val){
            this->real = val;
            this->imag = 0;
            return *this;
        }
                
        Complex& operator=(const int&val){
            this->real = val;
            this->imag = 0;
            return *this;
        }

        // bool operator<(const Complex&rhs){
        //     return fabs(imag) < TOL; 
        // }

        friend ostream&operator<<(ostream&out, Complex&obj){
            out<<setiosflags(ios::showpos|ios::fixed)<<setprecision(8)<<obj.real<<obj.imag<<"i";
            return out;
        }
};

//计算二阶矩阵的特征值
template<typename T1>
void root(T1 a,T1 b,T1 c,T1 d,Complex &s1, Complex&s2){
    T1 det = a*d - b*c;
    double delta = a*a + d*d + 2*a*d - 4*det;
    if(delta >=0){
        delta = sqrt(delta);
        s1 = Complex((a+d - delta)*0.5,0);
        s2 = Complex((a+d + delta)*0.5,0);
    }else{
        delta = sqrt(-delta);
        s1 = Complex((a+d)*0.5, -delta*0.5);
        s2 = Complex((a+d)*0.5, delta*0.5);
    }
}


//矩阵和向量模板类
template <typename T>
class Matrix{
    private:
        int row,col,size;
        T *pmm;
    public:
        Matrix(int r,int c):row(r),col(c){
            size = r*c;
            if(size >0){
                pmm = new T[size];
                for(int j=0;j<size;j++){
                    pmm[j] = 0;
                }
            } else
                pmm = NULL;
        }
        Matrix(int r,int c,T val):row(r),col(c){ //c为初始值
            size = r*c;
            if(size >0){
                pmm = new T[size];
                for(int j=0;j<size;j++){
                    pmm[j] = val;
                }
            } else
                pmm = NULL;
        }
        Matrix(const Matrix&rhs){ //copy constructor
            row = rhs.row;
            col = rhs.col;
            size = rhs.size;
            pmm = new T[size];
            for(int i=0;i<size;i++)
                pmm[i] = rhs.pmm[i];
        }
        ~Matrix(){
            if(!pmm){
                delete []pmm;
                pmm = NULL;
            }
        }
        Matrix<T>& operator=(const Matrix<T>&);
        int Row()const{return row;}
        int Col()const{return col;}
        friend ostream &operator<< <>(ostream&, const Matrix<T>&);          // 输出到屏幕
        friend Matrix<T> operator* <>(const Matrix<T>&lm,  const Matrix<T>&rm);
        friend Matrix<T> operator* <>(const Matrix<T>&lm, T val);
        friend Matrix<T> operator/ <>(const Matrix<T>&, T);
        friend Matrix<T> operator- <>(const Matrix<T>&lm,  const Matrix<T>&rm);
        friend Matrix<T> operator+ <>(const Matrix<T>&lm,  const Matrix<T>&rm);
        T*operator[](int i) const{return pmm + i*col;}
        Matrix<T> eye();
        Matrix<T> Transpose()const;
        void resize(int m,int n);
        // void sort();//排序
        Matrix<T> Hessen();  //变为上三角阵
        Matrix<Complex> DSQR(); //双步位移QR（double shift QR）
        void DQR(Matrix<T> &B); //双步位移使用的QR分解
        void QR(Matrix<T> &Q,Matrix<T> &R); //QR分解
        Matrix<T> eigVector(T v); //求特征值为v的特征向量
        double norm2();  //计算2范数
};



template<typename T>
double Matrix<T>::norm2(){
    double norm2 = 0;
    for(int i=0;i<size;i++){
        norm2 += pmm[i]*pmm[i];
    }
    return sqrt(norm2);
}

//生成对角阵
template<typename T>
Matrix<T> Matrix<T>::eye(){
    for(int i=0;i<row;i++){
        pmm[i*col +i] = 1;
    }
    return *this;
}

//转置
template<typename T>
Matrix<T> Matrix<T>::Transpose()const{
    Matrix<T> temp(col,row);
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
            temp[j][i] = pmm[i*col + j];
        }
    }
    return temp;
}

//赋值
template<typename T>
Matrix<T>& Matrix<T>::operator = (const Matrix<T>&rhs){
    if(this != &rhs){
        row = rhs.row;
        col = rhs.col;
        size = rhs.size;
        if(!pmm)
            delete []pmm;
        pmm = new T[size];
        for(int i=0;i<size;i++)
            pmm[i] = rhs.pmm[i];
    }
    return *this;
}

//降阶
template<typename T>
void Matrix<T>::resize(int m,int n){
    assert(row >=m && col >= n);
    for(int i=0;i<size;i++){
        pmm[i] = pmm[i/n*row + i%n];
    }
    row = m;col = n;
    size = m*n;
}

//变为Hessen矩阵，原始矩阵会改变
template<typename T>
Matrix<T> Matrix<T>::Hessen(){
    assert(row==col);
    double dr,cr,hr,tr;
    Matrix<T> a = *this;
    int n = a.row;
    Matrix<T> ur(n,1),pr(n,1),qr(n,1),wr(n,1); //看做只有一列的矩阵

    for(int r=0;r<n-2;r++){
        //1全为0的话Ar+1 = Ar
        bool flag = false;
        for(int i=r+2;i<n && !flag;i++){
            flag = fabs(a[i][r]) > TOL;
        }
        if(!flag) continue;
        //2计算一些常数
        double dr2 = 0;
        for(int i=r+1;i<n;i++){
            dr2 += a[i][r]*a[i][r];
        }
        dr = sqrt(dr2);
        cr = (a[r+1][r] > 0 ? -1 : 1)*dr;
        hr = dr2 -cr*a[r+1][r];
        //3计算ur
        for(int i=0;i<n;i++){
            if(i<r+1){
                ur[i][0] = 0;
            }else if(i== r+1){
                ur[i][0] = a[i][r] - cr;
            }else{
                ur[i][0] = a[i][r];
            }
        }
        //4计算Ar+1
        pr = (a.Transpose())*ur/hr;
        qr = a*ur*(1/hr);
        tr = (pr.Transpose()*ur/hr)[0][0];
        wr = qr - ur*tr;
        a = a - wr*(ur.Transpose()) - ur*(pr.Transpose());
    }
    cout<<"  Q2>对A拟上三角化之后的矩阵A(n-1)为："<<endl<<a<<endl;
    return a;
}

//QR分解
//输入参数为Q和R
template<typename T>
void Matrix<T>::QR(Matrix<T> &Q, Matrix<T> &A){  //此处A即为得到的R
    int m = A.row;
    Q = Q.eye();A = *this;
    double dr,cr,hr;
    Matrix<T> ur(m,1),pr(m,1),wr(m,1);
    for(int r=0;r<m-1;r++){
        int i=r+1;
        //1
        while(i<m && fabs(A[i][r])<= TOL)
            i++;
        if(i == m) continue;
        //2
        double dr2 = 0;
        for(int i=r;i<m;i++){
            dr2 += A[i][r]*A[i][r];
        }
        dr = sqrt(dr2);
        cr = (A[r][r] > 0 ? -1 : 1)*dr;
        hr = dr2 -cr*A[r][r];
        //3
        for(int i=0;i<m;i++){
            ur[i][0] = i<r ? 0 : i == r ? A[r][r] -cr : A[i][r];
        }
        //4
        wr = A*ur;
        Q = Q - wr*ur.Transpose()*(1/hr);
        pr = A.Transpose()*ur*(1/hr);
        A = A - ur*pr.Transpose();
    }
}

//双步位移中使用的QR分解
template<typename T>
void Matrix<T>::DQR(Matrix<T> &B){
    Matrix<T> &C = *this;
    int m = C.row;
    double dr,cr,hr,tr;
    Matrix<T> ur(m,1),vr(m,1),pr(m,1),qr(m,1),wr(m,1);
    for(int r=0;r<m-1;r++){
        int i=r+1;
        //1
        while(i<m && fabs(B[i][r])<= TOL) i++;
        if(i == m) continue;
        //2
        double dr2 = 0;
        for(int i=r;i<m;i++){
            dr2 += B[i][r]*B[i][r];
        }
        dr = sqrt(dr2);
        cr = (B[r][r] > 0 ? -1 : 1)*dr;
        hr = dr2 -cr*B[r][r];
        //3
        for(int i=0;i<m;i++){
            ur[i][0] = i<r ? 0 : i == r ? B[r][r] -cr : B[i][r];
        }
        //4
        vr = B.Transpose()*ur/hr;
        B = B - ur*vr.Transpose();
        pr = C.Transpose()*ur/hr;
        qr = C*ur/hr;
        tr = (pr.Transpose()*ur/hr)[0][0];
        wr = qr - ur*tr;
        C = C - wr*ur.Transpose() - ur*pr.Transpose();
    }
}

//双步位移计算特征值
template<typename T>
Matrix<Complex> Matrix<T>::DSQR(){
    Matrix<T> a = this->Hessen();
    int m = a.row;

    Matrix<T> Mk(m,m),I(m,m);
    Matrix<Complex> eig(m,1);
    I = I.eye();
    T s,t;

    int step = 0;
    while(m>2 && step <MAX_STEPS){
        //3
        if(fabs(a[m-1][m-2]) <= TOL){
            eig[m-1][0] = a[m-1][m-1];
            // cout<<"降1阶"<<a[m-1][m-1]<<" "<<eig[m-1][0]<<endl;
            m--;
            a.resize(m,m);I.resize(m,m);Mk.resize(m,m);
            continue;
        }     
        //
        if(fabs(a[m-2][m-3]) <= TOL){
            if(m == row){
                cout<<"  Q3.1>A(n-1)进行QR分解后得到的矩阵为:"<<endl;
                cout<<a<<endl;
                cout<<"  Q3.2>A(n-1)进行QR分解后的RQ值为："<<endl;
                cout<<(a*Mk)<<endl;
            }

           root(a[m-2][m-2],a[m-2][m-1],a[m-1][m-2],a[m-1][m-1],eig[m-2][0],eig[m-1][0]);
        //    cout<<"降2阶"<<eig[m-2][0]<<"  "<<eig[m-1][0]<<endl;
           m -=2;
           a.resize(m,m);I.resize(m,m);Mk.resize(m,m);
           continue; 
        }
        //9
        s = a[m-2][m-2] + a[m-1][m-1];
        t = a[m-2][m-2]*a[m-1][m-1] - a[m-1][m-2]*a[m-2][m-1];
        Mk = a*a - a*s + I*t;
        a.DQR(Mk);
        step++;
    }
    //4
    if(m == 1){
        eig[m-1][0] = a[m-1][m-1];
    }
    //5,6===
    if(m == 2){
        root(a[m-2][m-2],a[m-2][m-1],a[m-1][m-2],a[m-1][m-1],eig[m-2][0],eig[m-1][0]);
    }
    cout<<"  Q4>矩阵的全部特征值为:"<<endl;
    cout<<eig<<endl;
    return eig;
}

//计算特征向量
template<typename T>
Matrix<T> Matrix<T>::eigVector(T v){
    int n = row,index=n-1;
    Matrix<T> B(n,n), I(n,n),Q(n,n),R(n,n),V(n,1);
    I = I.eye();
    B = *this -I*v;
    B.QR(Q,R);
    //寻找特征值为0的位置
    for(int i=0;i<n;i++){
        if(fabs(R[i][i])<=TOL){
            index = i;
            break;
        }
    }
    //大于index的vi为0
    for(int i=index+1;i<n;i++){
        V[i][0] = 0;
    }
    //等于index的特征值先设置为1
    V[index][0] = 1;
    for(int i=index-1;i>=0;i--){
        for(int j=i+1;j<=index;j++){
            V[i][0] += R[i][j]*V[j][0];
        }
        V[i][0] /= -R[i][i];
    }
    V = V/V.norm2();
    return V;
}

#endif