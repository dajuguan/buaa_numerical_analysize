//hessen形式的qr分解法
#include<iostream>
#include<assert.h>
#include<math.h>
#include "matrix.h"
using namespace std;

#define N 10

void test_hessen(){
    Mat m(3,3);
    m[0][0] = 2;m[0][1] = 3;m[0][2] = 2;
    m[1][0] = 10;m[1][1] = 3;m[1][2] =4;
    m[2][0] = 3;m[2][1] = 6;m[2][2] =1;
    cout<<m<<endl;
    m.Hessen();
    cout<<"res:"<<m<<endl;
}

void test_root(){
    double a = 1,b=2,c=3,d=4;
    Complex s1,s2;
    double e = 0.00001;
    root(a,b,c,d,s1,s2);
    assert(fabs(s1.real + 0.372281)<e && fabs(s2.real - 5.37228)< e);
}

void test_DSQR(){
    Mat m(3,3);
    Matrix<Complex> res(3,3);
    m[0][0] = 2;m[0][1] = 3;m[0][2] = 2;
    m[1][0] = 10;m[1][1] = 3;m[1][2] =4;
    m[2][0] = 3;m[2][1] = 6;m[2][2] =1;
    res = m.DSQR();
    double e = 0.00001;
    assert(fabs(res[0][0].real - 11)<e && fabs(res[1][0].real +3)<e && fabs(res[2][0].real +2)<e);
}

//初始化矩阵
void genMat(Mat &m){
    int row = m.Row(),col = m.Col();
    for(int ii=0,i=ii+1;ii<row;ii++,i++){
        for(int jj=0,j=jj+1;jj<col;j++,jj++){
            m[ii][jj] = i==j ? 1.52*cos(i+1.2*j) : sin(0.5*i + 0.2*j);
        }
    }
}

//测试QR函数
void test_QR(){
    Mat m(3,3),q(3,3),r(3,3);
    m[0][0] = 2;m[0][1] = 3;m[0][2] = 2;
    m[1][0] = 10;m[1][1] = 3;m[1][2] =4;
    m[2][0] = 3;m[2][1] = 6;m[2][2] =1;   
    m.QR(q,r);
}
//解题程序
void solver(){
    Mat a(N,N),q(N,N),r(N,N);
    Matrix<Complex> eig(N,1);
    genMat(a);
    eig = a.DSQR();
    Mat v(N,1);
    for(int i=0;i<N;i++){
        if(fabs(eig[i][0].imag) < TOL){
            v = a.eigVector(eig[i][0].real); 
            cout<<"  Q5>实特征值为"<<eig[i][0].real<<"的特征向量为:"<<endl;
            cout<<v.Transpose()<<endl;    
        }
    }
}

int main(int argc, char const *argv[])
{
    solver();
}
