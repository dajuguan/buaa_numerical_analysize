//数值分析大作业第三题

//Copyright (c) 2018 chenjunfeng <wwwcjf163@163.com>
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.
//源码在https://github.com/dajuguan/buaa_numerical_analysize/tree/master/ex3

#include <iostream>
#include "matrix.h"
#include <math.h>

#define M 11
#define N 21 
#define TOL 1e-12
#define SIGMA 1e-7
#define MAX_STEPS 300

using namespace std;

const double mat_zuv[6][6]={
	{-0.5, -0.34, 0.14, 0.94, 2.06, 3.5},
	{-0.42, -0.5, -0.26, 0.3, 1.18, 2.38},
	{-0.18, -0.5, -0.5, -0.18, 0.46, 1.42},
	{0.22, -0.34, -0.58, -0.5, -0.1, 0.62},
	{0.78, -0.02, -0.5, -0.66, -0.5, -0.02},
	{1.5, 0.46, -0.26, -0.66, -0.74, -0.5}
};

const double arr_t[6] = {0, 0.2, 0.4, 0.6, 0.8, 1.0};
const double arr_u[6] = {0, 0.4, 0.8, 1.2, 1.6, 2.0};

//交换元素
void swap(double &a, double &b){
    double temp = a;
    a = b;
    b = temp;
}


/*选主元的杜立特尔分解法求线性方程组的解向量<组>
 *A为方程组系数矩阵
 *X为求解的方程组解向量组，按列排列（一个向量为1列）
 *B为方程组右端向量组，按列排列 
 * 
*/
void LU_Linea(Mat &A, Mat &X,Mat &B){
    assert( A.Col() == X.Row() && X.Row() == B.Row() && X.Col() == B.Col());
    int n = A.Col(),m = X.Col();
    Mat L(n,n);
    double s[n],maxSi;
    double y[n];
    int maxIndex;
    //1分解QA=LU
    for(int k=0;k<n;k++){
        maxSi = 0;maxIndex = 0;
        //计算中间量

        for(int i=k;i<n;i++){
            s[i] = A[i][k];
            for(int t=0;t<k;t++){
                s[i] -= L[i][t]*L[t][k];
            }
            //选取最大行号
            if(maxSi < fabs(s[i])){
                maxSi = fabs(s[i]);maxIndex = i;
            } 
        }
        if(maxIndex!=k){
            //交换
            L.swap(k,maxIndex,0,k);
            A.swap(k,maxIndex,k,n);
            B.swap(k,maxIndex,0,m);
            swap(s[k],s[maxIndex]);
        }
        //计算
        L[k][k] = s[k];
        for(int j=k+1;j<n&&k<n-1;j++){
            L[k][j] = A[k][j];
            for(int t=0;t<k;t++){
                L[k][j] -= L[k][t]*L[t][j];
            }
            L[j][k] = s[j]/L[k][k];
        }
    }        
    //3求解LY = Qb和Ux = y
    for(int r=0;r<m;r++){
        for(int i=0;i<n;i++){
            y[i] = B[i][r];
            for(int t=0;t<i;t++){
                y[i] -= L[i][t]*y[t];
            }
        }
        for(int i=n-1;i>=0;i--){
            X[i][r] = y[i];
            for(int t=i+1;t<n;t++){
                X[i][r] -= L[i][t]*X[t][r];
            }
            X[i][r] /= L[i][i];
        }
    }
}


//A为输入原矩阵，X为输出的逆矩阵
Mat inv(Mat A){
    int r = A.Row(), c = A.Col();
    Mat X(r,r), B(r,r);
    B.eye();
    LU_Linea(A,X,B);
    return X;
}


//生成题中的函数f和导函数ff
void genF(Vector &f,Mat &ff,double u,double v,double w,double t,double x,double y){
        f[0]=-1.0*(0.5*cos(t)+u+v+w-x-2.67);
        f[1]=-1.0*(t+0.5*sin(u)+v+w-y-1.07);
        f[2]=-1.0*(0.5*t+u+cos(v)+w-x-3.74);
        f[3]=-1.0*(t+0.5*u+v+sin(w)-y-0.79);
        ff[0][0]=-0.5*sin(t);
        ff[0][1]=1.0;
        ff[0][2]=1.0;
        ff[0][3]=1.0;
        ff[1][0]=1.0;
        ff[1][1]=0.5*cos(u);
        ff[1][2]=1.0;
        ff[1][3]=1.0;
        ff[2][0]=0.5;
        ff[2][1]=1.0;
        ff[2][2]=-sin(v);
        ff[2][3]=1.0;
        ff[3][0]=1.0;
        ff[3][1]=0.5;
        ff[3][2]=1.0;
        ff[3][3]=cos(w);
}

//牛顿法求解非线性方程组
void newton(const double x,const double y, double &t, double &u){
    Mat ff(4,4);
    Vector f(4),deltaX(4);
    double v,w;
    t = u = v = w = 1.0;
    int step = 0;
    while(step==0 ||(deltaX.absMax()/max(fabs(u),max(fabs(v),max(fabs(w),fabs(t)))) > TOL && step < MAX_STEPS)){
        genF(f,ff,u,v,w,t,x,y);
        LU_Linea(ff,deltaX,f);
        t+=deltaX[0];
        u+=deltaX[1];
        v+=deltaX[2];
        w+=deltaX[3];
        step++;
    }
    if(step == MAX_STEPS){
        cout<<"牛顿迭代法计算不成功"<<endl;
        exit;
    }
}

//分片双二次插值
double interp22(double x,double y){
	int i,j;
	double z = 0.0;
    i=int(fabs((x/0.2)+0.5));
    j=int(fabs((y/0.4)+0.5));
	if(i==0)  i=1;
    if(i==5)  i=4;
    if(j==0)  j=1;
    if(j==5)  j=4;

	for(int k=i-1;k<=i+1; k++)
		for(int r=j-1;r<=j+1;r++)
		{
			double sum=1.0;
			sum *= mat_zuv[k][r];
			for(int t1=i-1;t1<=i+1;t1++)
				if(t1!=k)
					sum*=(x-arr_t[t1])/(arr_t[k]-arr_t[t1]);
			for(int t2=j-1;t2<=j+1;t2++)
				if(t2!=r)
					sum*=(y-arr_u[t2])/(arr_u[r]-arr_u[t2]);
			z += sum;
		}	
	return z;
}


//对曲面拟合
double curveFitting(int k,double *X,double *Y,Mat &F,Mat &Crs){
    int n = F.Row();
    double sigma=0.0,p;
    Mat B(M,k+1),G(N,k+1),C(k+1,k+1),BT(k+1,M),GT(k+1,N);
    for(int i=0;i<M;i++)
        for(int j=0;j<=k;j++)
            B[i][j] = pow(0.08*i,j);
	for(int i=0;i<N;i++)
		for(int j=0;j<=k;j++)
            G[i][j]=pow(0.5+0.05*i,j);
    BT = B.Transpose();GT = G.Transpose();
    C = inv(BT*B);
    C = inv(BT*B)*BT*F*G*inv(GT*G);
    for(int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            p = 0.0;
            for(int r=0;r<=k;r++)
                for(int s=0;s<=k;s++)
                    p +=C[r][s]*pow(X[i],r)*pow(Y[j],s);
            sigma += (p - F[i][j])*(p - F[i][j]);
        };
    }
    cout<<"k="<<k<<"\t"<<"σ="<<setiosflags(ios::scientific)<<setprecision(12)<<sigma<<resetiosflags(ios::scientific)<<endl;
    if(sigma < SIGMA){
        cout<<"Q4>达到精度要求时的系数矩阵Crs:"<<endl;
        cout<<C<<endl;
    }
    Crs = C;
    return sigma;
}

//根据SIGMA要求选择k
void optimizeCurveFitting(double *X,double *Y, Mat &F,Mat &Crs){
    int k=2;
    cout<<"Q3>k的选择过程的k和σ:"<<endl;
    while(curveFitting(k,X,Y,F,Crs) > SIGMA) k++;
}

//计算拟合曲面输出的值
double fit(double x,double y,Mat &Crs){
    int k = Crs.Row();
    double p=0.0;
    for(int r=0;r<k;r++)
        for(int s=0;s<k;s++)
            p += Crs[r][s]*pow(x,r)*pow(y,s);    
    return p;
}

//测试函数用例
void testSwap(){
    Mat A(3,3),L(3,3),U(3,3);
    A[0][0] = 8.1;A[0][1] = 2.3; A[0][2] = -1.5;
    A[1][0] = 0.5;A[1][1] = -6.23;A[1][2] = 0.87;
    A[2][0] = 2.5;A[2][1] = 1.5; A[2][2] = 10.2;
    // LU_Linea(A,L);
    A.swap(0,2,0,2);
    cout<<A;
}

//测试函数用例
void testLU_Linea(){
    Mat A(2,2), L(2,2);
    Vector b(2),x(2);
    A[0][0] = 1e-20;A[0][1] = 1;
    A[1][0] = 1; A[1][1] = 2;
    b[0] = 1; b[1] = 4;
    LU_Linea(A,b,x);
    for(int i=0;i<2;i++){
        cout<<x[i]<<" ";
    }
    cout<<endl;
}

//测试函数用例
void testInv(){
    Mat A(3,3),X(3,3);
    A[0][0] = 8.1;A[0][1] = 2.3; A[0][2] = -1.5;
    A[1][0] = 0.5;A[1][1] = -6.23;A[1][2] = 0.87;
    A[2][0] = 2.5;A[2][1] = 1.5; A[2][2] = 10.2;
    X = inv(A);
    cout<<X<<endl;
}

//主程序
void solveF(){
    double u,t;
    Mat F(M,N),Crs(0,0);
    double x[M],y[N];
	for(int i=0;i<M;i++)
		x[i]=0.08*i;
	for(int j=0;j<N;j++)
		y[j]=0.5+0.05*j;
    //生成f(x,y)
    cout<<"Q2>数表x,y,f(x,y)为:"<<endl;
    for(int i=0;i<M;i++){
        for(int j=0;j<N;j++){
            newton(x[i],y[j],t,u);
            F[i][j] = interp22(t,u);
            cout<<noshowpos<<"x"<<i<<"="<<x[i]<<"\t"<<"  y"<<j<<"="<<fixed<<setprecision(2)<<y[j]<<"\t"
            <<"  f="<<resetiosflags(ios::fixed)<<setiosflags(ios::showpos|ios::scientific)<<setprecision(12)<<F[i][j]
			<<resetiosflags(ios::scientific)<<endl;
        }
    }
    //曲面拟合
    optimizeCurveFitting(x,y,F,Crs);
    //计算数表
    cout<<"Q5>数表x,y,f(x,y),p(x,y)为:"<<endl;
    for(int i=1;i<9;i++){
        for(int j=1;j<6;j++){
            double xx = 0.1*i, yy= 0.5+0.2*j,f,p;
            newton(xx,yy,t,u);
			cout<<std::noshowpos<<"x"<<i<<"*= "<<xx<<"\t"
			<<"y"<<j<<"*= "<<yy<<"\t"<<scientific<<showpos<<setprecision(12)<<"f="<<interp22(t,u)
			<<"\t"<<"p="<<fit(xx,yy,Crs)<<resetiosflags(ios::scientific)<<endl;  
        }
    }

}


int main(){
    solveF();
    return 1;
}