## Median filter by GuroMulgae ##

# Only 2D gray scale images

# %matplotlib notebook
# import numpy as np
# import cv2
# import matplotlib
# from matplotlib import pyplot as plt

def medfilt(A, n):

    # A : input image
    # n : n x n filter size
    
    A_mf = A.copy()
    
    nn = (n-1)//2
    [row, col] = CT_img.shape 
    ar = np.zeros((n,n))
    for r in range(row):
        for c in range(col):
            for j in range(-nn, nn + 1):
                for k in range(-nn, nn + 1):
                    if r+j <= -nn or r+j >= row + nn -1 or c+j <= -nn or c+j >= row + nn -1 or r+k <= -nn or r+k >= row + nn -1 or c+k <= -nn or c+k >= row + nn -1:
                        ar[j+1][k+1] = 0
                    else:
                        ar[j+1][k+1] = A[r+j][c+k]
            
            arr = np.sort(np.array(ar).ravel())
            A_mf[r][c] = arr[len(arr)//2]
    
    A_cv2_mf = cv2.medianBlur(A, 3)
    
    plt.figure(), plt.imshow(A, cmap='gray')
    plt.figure(), plt.imshow(A_mf, cmap='gray')
    plt.figure(), plt.imshow(A_cv2_mf, cmap='gray')
    
    
