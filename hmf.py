# GuroMulgae

def enh_hybridMedian(im,n=5):

    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    # Derive indices for the two patterns representing X and +
    indicesC = [0,4,6,8,12,16,18,20,24]
    indicesP = [2,7,10,11,12,13,14,17,22]
    
    v = (n-1) / 2

	# Process the image (ignoring the outer two layers of the image boundary
    for i in range(2,im.shape[0]-2):
        for j in range(2,im.shape[1]-2):
            # Extract the neighbourhood area
            block = im[i-v:i+v+1, j-v:j+v+1]
            
            # Reshape the neighborhood into a vector by flattening the 2D block
            wB = block.flatten()
            
            # Extract pixel values using indices
            wBc = numpy.take(wB,indicesC)
            wBp = numpy.take(wB,indicesP)
                  
            # Calculate the median values      
            wBcMed = numpy.median(wBc)
            wBpMed = numpy.median(wBp)
            
            # Calculate the hybrid median of the original pixel, and the two 
            # medians extracted above
            xmed = numpy.median([wBcMed,wBpMed,im[i][j]])

            # Assign the values               
            if (xmed > 0):
                img[i][j] = int(xmed)
            else:
                img[i][j] = im[i][j]
    return img
