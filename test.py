import numpy as np
#create matrix L with 4 rows and 6 columns with data:
#1 2 3 4 5 6
#7 8 9 10 11 12
#13 14 15 16 17 18
#19 20 21 22 23 24

L = np.array([[0,0,0,0,1,0],
              [0,0,0,0,0,1],
              [1,0,0,0,0,0],
              [0,1,0,0,0,0]])

K_CA = np.array([[1,2,3,4],
                [5,6,7,8],
                [9,10,11,12],
                [13,14,15,16]])

K = L.T @ K_CA @ L

print(K)