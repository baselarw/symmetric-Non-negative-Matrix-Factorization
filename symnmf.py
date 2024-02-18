import math
import sys
import numpy as np
import mysymnmf as smf


def symnf_func(k, goal, file):
    n=0     ## n= X.length()
    d=0     ## d= X[0].length()
    X =[]   ## X is an array of all the points
    my_file = open(file, "r")
    for line in my_file:
        if line != "\n":
            curr_victor = line.split(",")
            n+=1
            floating = []
            for j in curr_victor:
                floating.append(float(j))
            X.append(floating.copy())
    my_file.close()

    d = len(X[0])
    if(goal == 'sym'):
        res = smf.sym(X,n,d)
    elif(goal == 'ddg'):
        res = smf.ddg(X,n,d)
    elif(goal == 'norm'):
        res = smf.norm(X,n,d)
    elif(goal == 'symnmf'):
        W = smf.norm(X,n,d)
        H = initialize_H(W,n,k)
        res = smf.symnmf(H,W,n,k)

    ### printing the final result
    print_matrix(res,n)




def print_matrix(mat,n):
    result = [[0 for j in range(len(mat[0]))] for i in range(n)]
    for i in range(n):
        for j in range(len(mat[0])):
            result[i][j] = "%.4f" % mat[i][j]

    for i in range(n):
        str = ""
        for j in range(len(mat[0])):
            str += result[i][j] +","
        print(str[:-1])


def initialize_H(W,n,k):
    m = 0
    np.random.seed(0)
    res = [[0 for j in range(k)] for i in range(n)]
    for i in range(len(W)):
        for j in range(len(W)):
            m += W[i][j]
    m = m/(n*n)

    for i in range(n):
        for j in range(k):
            res[i][j] = np.random.uniform(0,2*math.sqrt(m/k),1)

    return res



inputs = sys.argv
inputs_length = len(inputs)    ### [not important, K, goal, File]
if(inputs_length != 4):
    print("An Error Has Occurred")
else:
    symnf_func(int(inputs[1]),inputs[2],inputs[3])