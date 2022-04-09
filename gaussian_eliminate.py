import numpy as np


def gaussian_elimiate(matrix):
    m = len(matrix)
    n = len(matrix[0])
    print(m)
    print(n)
    # if there are less rows than columns, then there wouldn't be a linear dependence
    if m < n:
        print("error")
        raise Exception("not enough data") 
    # set a m length array to locate the rows with pivots
    pivot = [0]*m
   
    pivot_dict = {}
    
    for j in range(n):

        # looks for the pivot in the column
        for i in range(m):
            # if a 1 is found at the i,j value
            # print(i, j, matrix[i][j])
            if(matrix[i][j] == 1):
                pivot[i] = 1
                # records the corresponding value for j where the pivot occurs
                pivot_dict[j]=i
                # adds the two rows using mod 2 addition
                for k in range(0,j):
                    if (matrix[i][k] == 1):
                        for row in range(m):
                            matrix[row][k] = (matrix[row][j] + matrix[row][k])%2

                for k in range(j+1,n):
                    if (matrix[i][k] == 1):
                        for row in range(m):
                            matrix[row][k] = (matrix[row][j] + matrix[row][k])%2
                break
    print(matrix)
    # stores which rows are dependent
    
    ret_all = []
    for i in range(m):
        ret = []
        # if there is a row that isnt a pivot, then it finds all the 1's in the row
        #  and the corresponding rows to those ones
        if pivot[i] == 0:
            ret.append(i)
            for key, col_value in enumerate(matrix[i]):
                if col_value == 1:
                    ret.append(pivot_dict[key])
            ret_all.append(ret)
                        
    return ret_all

if __name__ == '__main__':
    m = 12
    n = 10
    
    
    matrix = [[np.random.choice([0, 1], p=[.8, .2]) for _ in range(n)] for _ in range(m)]

    # matrix = [[1, 0, 0, 1, 0, 1],[1, 1, 0, 0, 0, 1],[0, 1, 0, 0, 1, 1], \
    #     [1, 1, 1, 0, 0, 1], [1, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 1, 1]]

    print(matrix)
    print(gaussian_elimiate(matrix))
