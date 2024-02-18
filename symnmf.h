# ifndef SYMNMF_H_
# define SYMNMF_H_

/*
*this function called when the goal == 'sym'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the sym mat
*/
double** sym(double** X,int n,int d);

/*
*this function called when the goal == 'ddg'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the ddg mat
*/
double** ddg(double** X,int n,int d);

/*
*this function called when the goal == 'norm'
*@param X: matrix[in this case the Set of the all points]
*@param n: num of rows
*@param d: num of cols
*return: the norm mat
*/
double** norm(double** X,int n,int d);

/*
*this function called when the goal == 'symnmf'
*@param H: matrix[in this case the initialized H]
*@param W: matrix[in this case the matrix that we get from norm()]
*@param n: num of rows
*@param d: num of cols
*return: the final mat of the completed algro
*/
double** symnmf(double** H,double** W,int n,int k);

/*
*this function called to free the allocated memory [in this project all of them are matrices]
*@param mat: matrix
*@param rows: num of rows
*/
void freeMatrix(double **matrix, int rows);
# endif
