def gaussian_elimiate(matrix):
    m = len(matrix)
    n = len(matrix[0])
    # if there are less rows than columns, then there wouldn't be a linear dependence
    if(m < n):
        raise Exception("not enough data") 
    # set a m length array to locate the rows with pivots
    pivot = [0]*m
   
    pivot_dict = {}
    
    for j in range(n):

        # looks for the pivot in the column
        for i in range(m):
            # if a 1 is found at the i,j value
            if(matrix[i][j] == 1):
                pivot[i] = 1
                # records the corresponding value for j where the pivot occurs
                pivot_dict[j]=i
                # adds the two rows using mod 2 addition
                for k in range(0,i):
                    if (m[i][k] == 1):
                        for row in range(m):
                            m[row][k] = (m[row][j] + m[row][k])%2

                for k in range(i+1,m):
                    if (m[i][k] == 1):
                        for row in range(m):
                            m[row][k] = (m[row][j] + m[row][k])%2
            break
    
    # stores which rows are dependent
    ret = []
    for i in range(m):
        # if there is a row that isnt a pivot, then it finds all the 1's in the row
        #  and the corresponding rows to those ones
        if pivot[i] == 0:
            ret.append[i]
            for key, col_value in enumerate(matrix[i]):
                if col_value == 1:
                    ret.append(pivot_dict[key])
                        
    return ret