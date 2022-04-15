import numpy as np

def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)


def gaussian_elimiate(matrix2):
    matrix_copy = matrix2
    m = len(matrix2)
    n = len(matrix2[0])
    
    # transposes for faster calculation
    matrix = transpose(matrix2)

    # print(matrix)
    # print(matrix)
    # if there are less rows than columns, then there wouldn't be a linear dependence
    if m < n:
        print("error")
        raise Exception("not enough data") 
    # set a m length array to locate the rows with pivots
    pivot = [0]*m
   
    pivot_dict = {}
    
    # iterates through all the columns of the original matrix (rows of new matrix)
    for i in range(n):

        # looks for the pivot in the column (in the row of the new matrix)
        for j in range(m):
            # # if a 1 is found at the i,j value
            # # print(i, j, matrix[i][j])
            # if(matrix[i][j] == 1):
            #     pivot[i] = 1
            #     # records the corresponding value for j where the pivot occurs
            #     pivot_dict[j]=i
            #     # adds the two rows using mod 2 addition
            #     for k in range(0,j):
            #         if (matrix[i][k] == 1):
            #             for row in range(m):
            #                 matrix[row][k] = (matrix[row][j] + matrix[row][k])%2

            #     for k in range(j+1,n):
            #         if (matrix[i][k] == 1):
            #             for row in range(m):
            #                 matrix[row][k] = (matrix[row][j] + matrix[row][k])%2
            #     break
            if(matrix[i][j] == 1):
                pivot[j] = 1
                # records the corresponding value for j where the pivot occurs
                pivot_dict[i]=j
                # adds the two rows using mod 2 addition
                for k in range(0,i):
                    if (matrix[k][j] == 1):
                        for row in range(m):
                            matrix[k][row] = (matrix[i][row] + matrix[k][row])%2

                for k in range(i+1,n):
                    if (matrix[k][j] == 1):
                        for row in range(m):
                            matrix[k][row] = (matrix[i][row] + matrix[k][row])%2
                break
    
    # stores which rows are dependent

    
    # print(matrix)
    # matrix = transpose(matrix)
    # print(matrix)
    # print(pivot_dict)
    # print(pivot)
    ret_all = []
    for j in range(m):
        ret = []
        # if there is a row that isnt a pivot, then it finds all the 1's in the row
        #  and the corresponding rows to those ones
        if pivot[j] == 0:
            ret.append(j)
            for i in range(len(matrix)):
                if matrix[i][j] == 1:
                    ret.append(pivot_dict[i])
            ret_all.append(ret)
                        
    return ret_all

if __name__ == '__main__':
    m = 6
    n = 5
    
    
    matrix = [[np.random.choice([0, 1], p=[.8, .2]) for _ in range(n)] for _ in range(m)]

    # matrix = [[1, 0, 0, 1, 0, 1],[1, 1, 0, 0, 0, 1],[0, 1, 0, 0, 1, 1], \
    #     [1, 1, 1, 0, 0, 1], [1, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 1, 1]]
    # matrix = [[0, 0, 0, 1, 1], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 1, 0]]
    print(matrix)
    print(gaussian_elimiate(matrix))
