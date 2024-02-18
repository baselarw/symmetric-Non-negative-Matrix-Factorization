#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

double** sym(double** X,int n,int d);
double eucldis(double *x1,double *x2,int d);
double** ddg(double** X,int n,int d);
double** norm(double** X,int n,int d);
double** symnmf(double** H,double** W,int n,int k);
double** mat_transpose(double** X,int rows, int cols);
void freeMatrix(double **matrix, int rows);
void printmatrix(double **mat,int rows,int cols);
double** multiplyMatrices(double** firstMatrix, double** secondMatrix, int rowFirst, int colFirst, int rowSecond, int colSecond);
int main(int argc, char **argv){
    int n=0,d=0,firstvec=1,i=0,j=0;
    double **X,**res;
    double temp;
    char c;
    FILE *firstscan=NULL,*secondscan=NULL;
    if(argc!=3){
        printf("An Error Has Occurred\n");
        return 1;
    }
    firstscan=fopen(argv[2],"r");
    if(firstscan==NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    while(fscanf(firstscan, "%lf%c", &temp,&c) == 2){
        if(firstvec==1){d+=1;}
        if(c!=','){
            n+=1;
            firstvec=0;
        }
    }
    fclose(firstscan);
    secondscan=fopen(argv[2],"r");
    if(secondscan==NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    X=(double **) malloc(n*sizeof(double *));
    for(i=0;i<n;i++){
        X[i]=(double *)malloc(d*sizeof(double));
    }
    i=0,j=0;
    while(fscanf(secondscan, "%lf%c", &temp,&c) >= 1){
        if(c!=EOF){
            if(c==','){
                X[i][j]=temp;
                j+=1;
            }
            else{
                X[i][j]=temp;
                i+=1;
                j=0;
            }
        }
        else{/*c=EOF*/
                X[i][j]=temp;
        }
    }
    if(strcmp(argv[1],"sym")==0){
        res=sym(X,n,d);
        printmatrix(res,n,n);
        freeMatrix(res,n);
    }
    if(strcmp(argv[1],"ddg")==0){
        res=ddg(X,n,d);
        printmatrix(res,n,n);
        freeMatrix(res,n);
    }
    if(strcmp(argv[1],"norm")==0){
        res=norm(X,n,d);
        printmatrix(res,n,n);
        freeMatrix(res,n);
    }
    freeMatrix(X,n);
    fclose(secondscan);
    return 0;
}

/*
*this function called when the goal == 'symnmf'
*@param H: matrix[in this case the initialized H]
*@param W: matrix[in this case the matrix that we get from norm()]
*@param n: num of rows
*@param d: num of cols
*return: the final mat of the completed algro
*/
double** symnmf(double** H,double** W,int n,int k){
    double eps=0.0001,beta=0.5,convergence=eps+3,temp=0;
    int max_iter=300,iter=0,i,j;
    double **HT,**WH,**HHTH;
    double **temp1;
    HT=mat_transpose(H,n,k);
    WH=multiplyMatrices(W,H,n,n,n,k);
    temp1=multiplyMatrices(H,HT,n,k,k,n);
    HHTH=multiplyMatrices(temp1,H,n,n,n,k);
    while(iter<max_iter && convergence>=eps){
        convergence=0;
        for(i=0;i<n;i++){
            for(j=0;j<k;j++){
                temp=H[i][j];
                H[i][j]=H[i][j]*(1-beta+((beta*WH[i][j])/HHTH[i][j]));
                convergence+=pow(temp-H[i][j],2);
            }
        }
        freeMatrix(temp1,n);
        freeMatrix(HT,k);
        freeMatrix(WH,n);
        freeMatrix(HHTH,n);
        HT=mat_transpose(H,n,k);
        WH=multiplyMatrices(W,H,n,n,n,k);
        temp1=multiplyMatrices(H,HT,n,k,k,n);
        HHTH=multiplyMatrices(temp1,H,n,n,n,k);
        iter+=1;
    }
    freeMatrix(HT,k);
    freeMatrix(WH,n);
    freeMatrix(HHTH,n);
    freeMatrix(temp1,n);
    return H;
}

/*
*@param X: matrix
*@param rows: num of rows
*@param cols: num of cols
*return: the transopse of the matrix
*/
double** mat_transpose(double** X,int rows, int cols){
    double **res=(double **) malloc(cols*sizeof(double *));
    int i,j;
    for(i=0;i<cols;i++){
        res[i]=(double *)malloc(rows*sizeof(double));
    }
    for(i=0;i<cols;i++){
        for(j=0;j<rows;j++){
            res[i][j]=X[j][i];
        }
    }
    return res;
}

/*
*this function called when the goal == 'sym'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the sym mat
*/
double** sym(double** X,int n,int d){
    double **res=(double **) malloc(n*sizeof(double *));
    int i,j;
    for(i=0;i<n;i++){
        res[i]=(double *)malloc(n*sizeof(double));
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){res[i][j]=0;}
            else{
                res[i][j]=exp(-0.5*eucldis(X[i],X[j],d));
            }
        }
    }
    return res;
}

/*
*this function calculate the eucldis distance
*@param x1: first point
*@param x2: second point
*@param d: num of cols (len(x1)==len(x2)==d)
*return: the eucldis distance
*/
double eucldis(double *x1,double *x2,int d){
    int i;
    double dis=0;
    for(i=0;i<d;i++){
        dis+= pow(x1[i]-x2[i],2);
    }
    return dis;
}

/*
*this function called when the goal == 'ddg'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the ddg mat
*/
double** ddg(double** X,int n,int d){
    double **res=(double **) malloc(n*sizeof(double *));
    double **symmat=sym(X,n,d);
    int i,j,col;
    double di;
    for(i=0;i<n;i++){
        res[i]=(double *)malloc(n*sizeof(double));
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            res[i][j]=0;
            if(i==j){
                di=0;
                for(col=0;col<n;col++){di+=symmat[i][col];}
                res[i][j]=di;
            }
        }
    }
    freeMatrix(symmat,n);
    return res;
}

/*
*this function called when the goal == 'norm'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the norm mat
*/
double** norm(double** X,int n,int d){
    double **A = sym(X,n,d);
    double **D = ddg(X,n,d);
    double **temp,**W;
    double **D_inv = (double **) malloc(n*sizeof(double *));
    int i,j;
    for(i=0; i<n; i++){
        D_inv[i] = (double *) malloc(n*sizeof(double));
    }

    for(i=0; i<n; i++){       
        for(j=0; j<n; j++){
            if(i==j){
                D_inv[i][j] = pow(D[i][j],-0.5);
            }
            else{
                D_inv[i][j] =0;
            }
        }
    }
    temp=multiplyMatrices(D_inv,A,n,n,n,n);
    W = multiplyMatrices(temp,D_inv,n,n,n,n);
    
    freeMatrix(temp,n);
    freeMatrix(A,n);
    freeMatrix(D,n);
    freeMatrix(D_inv,n);
    return W;
}

/*
*this function called to print the matrix
*@param mat: matrix
*@param rows: num of rows
*@param cols: num of cols
*/
void printmatrix(double **mat,int rows,int cols){
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            if(j!=cols-1){printf("%.4f,",mat[i][j]);}
            else{printf("%.4f\n",mat[i][j]);}
        }
    }
}

/*
*this function called to free the allocated memory [in this project all of them are matrices]
*@param mat: matrix
*@param rows: num of rows
*/
void freeMatrix(double **matrix, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/*
*this function called to print the matrix
*@param firstMatrix: matrix
*@param secondMatrix: matrix
*@param rowFirst: num of the rows of the first matrix
*@param colFirst: num of the cols of the first matrix
*@param rowSecond: num of the rows of the second matrix
*@param colSecond: num of the cols of the second matrix
*@return: the result of the multiply
*/
double** multiplyMatrices(double** firstMatrix, double** secondMatrix, int rowFirst, int colFirst, int rowSecond, int colSecond) {
    int i,j,k;
    double** result = (double**)malloc(rowFirst * sizeof(double*));
    if(firstMatrix==NULL || secondMatrix==NULL){
        freeMatrix(result,rowFirst);
        printf("An Error Has Occurred\n");
        return NULL;
    }
    if (colFirst != rowSecond) {
        freeMatrix(result,rowFirst);
        printf("An Error Has Occurred\n");
        return NULL;
    }
    for (i = 0; i < rowFirst; ++i) {
        result[i] = (double*)malloc(colSecond * sizeof(double));
    }
    for (i = 0; i < rowFirst; ++i) {
        for (j = 0; j < colSecond; ++j) {
            result[i][j] = 0.0;
            for (k = 0; k < colFirst; ++k) {
                result[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
            }
        }
    }
    return result;
}