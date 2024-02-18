import math
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import sys
import mysymnmf as smf


def symnf_func(k, file):
    n = 0  ## n= X.length()
    d = 0  ## d= X[0].length()
    X = []  ## X is an array of all the points
    my_file = open(file, "r")
    for line in my_file:
        if line != "\n":
            curr_victor = line.split(",")
            n += 1
            floating = []
            for j in curr_victor:
                floating.append(float(j))
            X.append(floating.copy())
    my_file.close()

    d = len(X[0])
    W = smf.norm(X, n, d)
    H = initialize_H(W, n, k)
    H_final_mat = smf.symnmf(H, W, n, k)
    ### calculating the nfm answer:
    nmf_arr = clusters_adjustment(H_final_mat, k)
    nmf_arr_np = np.array(nmf_arr)
    X_np = np.array(X)
    nmf_res = "%.4f" % silhouette_score(X_np,nmf_arr_np)

    ### calculating the kmeans answer:
    Kmeans_arr = K_means_func(k,300,X)
    Kmeans_arr_np = np.array(Kmeans_arr)
    Kmeans_res = "%.4f" % silhouette_score(X_np,Kmeans_arr_np)
    print("nmf: " + nmf_res)
    print("kmeans: " + Kmeans_res)


def K_means_func(k,iter,X):
    vectors_float = X
    means = []

    Kmeans_arr=[0 for i in range(len(X))]  # result to return

    epsilon = 0.0001
    for i in range(k):              ### adding the first k vectors as centroids///// and the i'th mean is fitting to the i'th cluster
        means.append(vectors_float[i].copy())
    l=0
    while ((l < iter)) :

        counters = [0 for i in range(k)]
        summ_vectors_cluster = [[] for i in range(k)]

        vectors_length = len(vectors_float)

        old_means = [[] for i in range(k)]    ### [[],[],[]]
        for i in range(k):
            for j in range(len(means[0])):
                old_means[i].append(means[i][j])

        for i in range (vectors_length):
            curr_vec = vectors_float[i].copy()
            min_dist = distance(curr_vec,means[0])
            min_dist_index = 0
            for m in range(1,k):             ## finding the min_dist_index
                curr_dist = distance(curr_vec, means[m])
                if(min_dist>curr_dist):
                    min_dist = curr_dist
                    min_dist_index = m
            Kmeans_arr[i] = min_dist_index   # updating the Kmeans_arr
            if(counters[min_dist_index] ==0):
                summ_vectors_cluster[min_dist_index] = curr_vec.copy()
            else:    ## adding to the old summ of the cluster another victor
                extract_vec = summ_vectors_cluster[min_dist_index]
                for n in range(len(curr_vec)):
                    extract_vec[n] += curr_vec[n]
            counters[min_dist_index] += 1
        ### now we will calc the new means:
        for i in range(k):
            for j in range(len(means[0])):
                means[i][j] =summ_vectors_cluster[i][j]/counters[i]

        ## checking the epsilon
        cnt=0
        for i in range(k):
            curr_dist = distance(old_means[i],means[i])
            if(curr_dist < 0.001):
                cnt+=1
        if(cnt == k):
            break
        l+=1
    return Kmeans_arr






def distance(victor1, victor2):         ### calcuting the Euclidean Distance
    n = len(victor1)
    result = 0
    for i in range(n):
        sub = victor1[i]-victor2[i]
        poww = pow(sub,2)
        result += poww

    return pow(result,0.5)

def clusters_adjustment(mat,k):
    arr_res = [0 for i in range(len(mat))]
    for i in range(len(mat)):
        max_point = mat[i][0]
        max_index = 0
        for j in range(1,k):
            curr_point = mat[i][j]
            if(curr_point>max_point):
                max_point = curr_point
                max_index = j
        arr_res[i]=max_index
    return arr_res




def initialize_H(W, n, k):
    m = 0
    np.random.seed(0)
    res = [[0 for j in range(k)] for i in range(n)]
    for i in range(len(W)):
        for j in range(len(W)):
            m += W[i][j]
    m = m / (n * n)

    for i in range(n):
        for j in range(k):
            res[i][j] = np.random.uniform(0, 2 * math.sqrt(m / k), 1)

    return res


inputs = sys.argv
inputs_length = len(inputs)  ### [not important, K, File]
if (inputs_length != 3):
    print("An Error Has Occurred")
else:
    symnf_func(int(inputs[1]), inputs[2])