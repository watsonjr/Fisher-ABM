function parallel_julia(I)  
   # determine what section of the matrix we're calculating  
   yrange = I[1]  
   xrange = I[2]  
   array_slice = (size(yrange, 1), size(xrange, 1))  

   # construct interim matrix that we'll populate  
   matrix = Array(Uint8, array_slice)  

   # determine the lowest coordinates  
   xmin=xrange[1]  
   ymin=yrange[1]  

   # get the pixel values, convert to uint8  
   for x=xrange, y=yrange  
     pix = 256
     z = complex((x-width/2)/(height/2), (y-height/2)/(height/2)) 

     # find the value of the pixel using the above algorithm 
     for n = 1:256
       if abs(z) > 2
         pix = n-1
         break
       end
       z = z^2 + C
     end
     matrix[y-ymin+1, x-xmin+1] = uint8(pix)
   end  

   # done!  
   return matrix  
 end  


