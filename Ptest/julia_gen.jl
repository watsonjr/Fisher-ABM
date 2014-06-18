# use the Images library  
 using Images, PyPlot
 # load up our matrix calculation function  
 require("julia_matrix.jl")  

 # define some variable to be available to ALL processors  
 @everywhere width = 1000  
 @everywhere height = 700  
 @everywhere C = -0.8 - 0.65im  

 # construct our distributed array, passing in our function  
 # and the dimensions of the full pixel matrix  
 distributed_matrix = DArray(parallel_julia, (height, width))  

 # convert it to an array we can write our output file  
 pixel_matrix = convert(Array, distributed_matrix)  
 pcolormesh(pixel_matrix)
 imwrite(pixel_matrix, "julia_set.png")  


