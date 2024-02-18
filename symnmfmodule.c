#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

double** PYmat_to_Cmat(PyObject *py_mat, int n, int d);
void Cmat_to_PYmat(double **c_mat,PyObject *py_mat, int rows, int colms);

static PyObject* sym_Wrapper(PyObject *self, PyObject *args)
{
    PyObject* X;
    int n,d;
    double** res_Cmat;
    double** py_C_mat;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "Oii",&X,&n,&d)) {
        return NULL;
    }
    py_C_mat = PYmat_to_Cmat(X,n,d);
    res_Cmat = sym(py_C_mat,n,d);
    res=PyList_New(n);
    Cmat_to_PYmat(res_Cmat,res,n,n);
    freeMatrix(res_Cmat,n);
    freeMatrix(py_C_mat,n);
    return res;
}

static PyObject* ddg_Wrapper(PyObject *self, PyObject *args)
{
    PyObject* X;
    int n,d;
    double** res_Cmat;
    double** py_C_mat;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "Oii",&X,&n,&d)) {
        return NULL;
    }
    py_C_mat = PYmat_to_Cmat(X,n,d);
    res_Cmat = ddg(py_C_mat,n,d);
    res=PyList_New(n);
    Cmat_to_PYmat(res_Cmat,res,n,n);
    freeMatrix(res_Cmat,n);
    freeMatrix(py_C_mat,n);
    return res;
}

static PyObject* norm_Wrapper(PyObject *self, PyObject *args)
{
    PyObject* X;
    int n,d;
    double** res_Cmat;
    double** py_C_mat;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "Oii",&X,&n,&d)) {
        return NULL;
    }
    py_C_mat = PYmat_to_Cmat(X,n,d);
    res_Cmat = norm(py_C_mat,n,d);
    res=PyList_New(n);
    Cmat_to_PYmat(res_Cmat,res,n,n);
    freeMatrix(res_Cmat,n);
    freeMatrix(py_C_mat,n);
    return res;
}

static PyObject* symnmf_Wrapper(PyObject *self, PyObject *args)
{
    PyObject *H;
    PyObject *W;
    int n,k;
    double** res_Cmat;
    double** py_C_mat_H;
    double** py_C_mat_W;
    PyObject* res;
    if(!PyArg_ParseTuple(args, "OOii",&H,&W,&n,&k)) {
        return NULL;
    }
    
    py_C_mat_H = PYmat_to_Cmat(H,n,k);
    py_C_mat_W = PYmat_to_Cmat(W,n,n);
    res_Cmat = symnmf(py_C_mat_H, py_C_mat_W, n, k);

    res=PyList_New(n);
    Cmat_to_PYmat(res_Cmat,res,n,k);
    freeMatrix(res_Cmat,n);
    freeMatrix(py_C_mat_W,n);
    return res;
}

static PyMethodDef SymNMFMethods[] = {
    {"sym",                   
      (PyCFunction) sym_Wrapper, 
      METH_VARARGS,      
      PyDoc_STR("the function gets the datapoints X, the number of datapoints n, and the dimension of each point and calculates the similarity matrix")
    }, 
    {"ddg",                   
      (PyCFunction) ddg_Wrapper, 
      METH_VARARGS,      
      PyDoc_STR("the function gets the datapoints X, the number of datapoints n, and the dimension of each point and calculates the diagonal degree matrix")
    },
    {"norm",                   
      (PyCFunction) norm_Wrapper, 
      METH_VARARGS,      
      PyDoc_STR("the function gets the datapoints X, the number of datapoints n, and the dimension of each point and calculates the normalized similarity matrix")
    },
    {"symnmf",                   
      (PyCFunction) symnmf_Wrapper, 
      METH_VARARGS,      
      PyDoc_STR("the function gets an initaial H, the W and two other arguments n,k when H:nxk and W:nxn anc calculates the final H")

    },
    {NULL, NULL, 0, NULL}   
};

static struct PyModuleDef Symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf", 
    NULL,
    -1,  
    SymNMFMethods 
};

PyMODINIT_FUNC PyInit_mysymnmf(void) {
    PyObject *m;
    m = PyModule_Create(&Symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

/*
*This function convert from python matrix to C matrix
*@param py_mat:python Object(in this case the python matrix)
*@param n: num of rows
*@param d: num of cols
*@return: C matrix (after converting from python matrix)
*/
double** PYmat_to_Cmat(PyObject *py_mat, int n, int d){
    int i,j;
    double **c_mat=(double **) malloc(n*sizeof(double *));
    PyObject *py_row;
    for (i=0;i<n;i++){
        c_mat[i] = (double*) malloc(d*sizeof(double));
        if(c_mat[i] == NULL){
            printf("An Error Has Occurred");
        }
    }
    for(i=0;i<n;i++){
        py_row = PyList_GetItem(py_mat, i);
        for (j=0;j<d;j++){
            c_mat[i][j] = PyFloat_AsDouble(PyList_GetItem(py_row, j));
        }
    }
    return c_mat;
}

/*
*This function convert from C matrix to python matrix
*@param py_mat:python Object(in this case the python matrix)
*@param C_mat:python Object(in this case the C matrix)
*@param n: num of rows
*@param d: num of cols
*Explanation: in the function above we return the matrix
* in this function we change the python matrix by the pointer
* two different ways to make the change
*/
void Cmat_to_PYmat(double **c_mat,PyObject *py_mat, int rows, int colms){
    PyObject *py_col_mat;
    int i,j;
    for(i=0;i<rows;i++){
        py_col_mat = PyList_New(colms);
        for(j=0;j<colms;j++)
        {
            PyList_SetItem(py_col_mat,j,PyFloat_FromDouble(c_mat[i][j]));
        }
        PyList_SetItem(py_mat,i,py_col_mat);
    }
}