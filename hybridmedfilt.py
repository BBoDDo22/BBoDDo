## Hybrid Median Filter by GuroMulgae ##

def hybridmedfilt(A, n):

    # A : input image
    # n : n x n kernel size
    
    # needs odd n
    n = round(n)
    if (n + 1) % 2 != 0:
        n = n + 1
        
    lineC1 = np.concatenate((np.arange(0, n**2-1+1, n+1), np.arange(n-1, n*(n-1)+1, n-1)))
    lineC = np.sort(np.delete(lineC1, np.where(lineC1 == (n**2-1)//2))) # two lines crossing to each other
    
    lineP1 = np.concatenate((np.arange((n-1)//2, n**2-(n-1)//2+1, n), np.arange(n*(n-1)//2, n*(n+1)//2-1+1, 1)))
    lineP = np.sort(np.delete(lineP1, np.where(lineP1 == (n**2-1)//2))) # two lines perpendicular to each other
    
    A_hmf = A.copy()
    
    nn = (n-1)//2
    [row, col] = A.shape 
    
    ar = np.zeros((n,n))
    for r in range(nn,row-nn):
        for c in range(nn,col-nn):
            # Extract the image nxn mask
            ar = A[r-nn:r+nn+1, c-nn:c+nn+1]
            
            # Flatten the 2D array to 1D
            arr = ar.ravel()
            
            # Extract line values 
            arr_c = np.take(arr, lineC)
            arr_p = np.take(arr, lineP)
            
            # Select median values
            arr_c_med = arr_c[len(arr_c)//2]
            arr_p_med = arr_p[len(arr_p)//2]
            
            # Pick median and result
            med = np.sort([arr_c_med, arr_p_med, A[r][c]])
            A_hmf[r][c] = med[len(med)//2]

    return A_hmf
