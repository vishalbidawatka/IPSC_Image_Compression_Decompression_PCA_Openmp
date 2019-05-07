# IPSC_Image_Compression_Decompression_PCA
## Team Number : 5
##### Vishal Bidawatka (2018201004)
##### Kshitij Paliwal (2018201063)
##### Sandeep Kumar Gupta (2018201076)

<br>
<br>
</br>


###### Link for downloading data :-
https://drive.google.com/open?id=1r8g20jgtmGH6lZ5gfF7ZoxLDrLTxok-s

### Steps to be followed for parallel code :

1. Download 0.txt or 8.txt from the given link.
2. Compile parallel_pca as : g++ -std=c++14 -fopenmp Parallel_pca.cpp -o svd 
3. run parallel_pca as ./svd [rows] [cols] [NUM_THREADS] [FILENAME_of_datamatrix] [Percentage_to_compress] 

### Example of datamatrix "8.txt" with 4 threads and 90% compression
- Rows and columns will be typically 784 as we have tested it on mnist dataset
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/screenshot.png)

### Steps to followed for serial code : 

1. Simply compile and run the serial_pca.


- The ouput after running any of these files is generated into reconstructed_matrix.txt.
- Open python preprocessing jupyter notebook and run the cell generating output from given file.

# Results
- Original image (0% compression )
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/exact_8.png)

- Reconstructed from 10% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/10_precent_8.png)

- Reconstructed from 20% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/20_precent_8.png)

- Reconstructed from 40% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/40_precent_8.png)

- Reconstructed from 60% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/60_precent_8.png)

- Reconstructed from 80% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/80_precent_8.png)

- Reconstructed from 90% compression
![Original Image ](https://github.com/vishalbidawatka/IPSC_Image_Compression_Decompression_PCA/blob/master/reresultimagezero/90_precent_8.png)
