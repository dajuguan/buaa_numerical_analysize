//本文件是数值分析第三个作业定义的通用矩阵和向量操作头文件matrix.h

#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>
#include<iomanip>
#include<assert.h>
#include<math.h>
using namespace std;

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
    out<<resetiosflags(ios::showpos|ios::fixed)<<endl;
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


//矩阵和向量模板类
template <typename T>
class Matrix{
    protected:
        T *pmm;
        int row,col,size;
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
        Matrix(int r,int c,T* arr):row(r),col(c){ //c为初始值
            size = r;
            if(size >0){
                pmm = new T[size];
                for(int j=0;j<size;j++){
                    pmm[j] = arr[j];
                }
            } else
                pmm = NULL;
        }
        Matrix(const Matrix&rhs){ //copy constructor
                row = rhs.row;
                col = rhs.col;
                size = rhs.size;
                if(size > 0){
                    pmm = new T[size];
                    for(int i=0;i<size;i++)
                        pmm[i] = rhs.pmm[i];
                } else pmm = NULL;

        }
        ~Matrix(){
            if(pmm){
                delete [] pmm;
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
        void eye();
        Matrix<T> Transpose()const;
        void resize(int m,int n);
        void swap(int i,int k,int start,int end); //交换两行元素
        // void sort();//排序
        Matrix<T> Hessen();  //变为上三角阵
        void DQR(Matrix<T> &B); //双步位移使用的QR分解
        void QR(Matrix<T> &Q,Matrix<T> &R); //QR分解
        Matrix<T> eigVector(T v); //求特征值为v的特征向量
        double norm2();  //计算2范数
        double absMax();  //绝对值最大的元素
};

class Vector:public Matrix<double>{
    public:
        Vector(int n):Matrix<double>(n,1){}
        double &operator[](int i) const{return *(pmm + i*col);}
};


template<typename T>
double Matrix<T>::norm2(){
    double res = 0;
    for(int i=0;i<size;i++){
        res += pmm[i]*pmm[i];
    }
    return sqrt(res);
}

//生成对角阵
template<typename T>
void Matrix<T>::eye(){
    for(int i=0;i<row;i++){
        pmm[i*col +i] = 1;
    }
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
        if(pmm) delete []pmm;
        pmm = new T[size];
        if(pmm){
            for(int i=0;i<size;i++)
                pmm[i] = rhs.pmm[i];
        }
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


//交换两行元素
template<typename T>
void Matrix<T>::swap(int m, int n,  int start, int end){
    assert( m>-1 && m <row && n>-1 && n < row && start > -1 && start<col  && end<= col);
    for(int i=start;i<end;i++){
        T temp = pmm[m*col + i];
        pmm[m*col + i] = pmm[n*col + i];
        pmm[n*col + i] = temp;
    }
}

//绝对值最大的元素
template<typename T>
double Matrix<T>::absMax(){
    double _max = 0.0;
    for(int i=0;i<size;i++)
        _max < fabs(pmm[i]) ? _max = fabs(pmm[i]) : NULL;
    return _max;
}

#endif