# IPSC_Image_Compression_Decompression_PCA

- Link for downloading data :-
https://drive.google.com/open?id=1r8g20jgtmGH6lZ5gfF7ZoxLDrLTxok-s

- Steps to be followed for parallel code :

1. Download 0.txt or 8.txt from the given link.
2. Compile parallel_pca as : g++ -std=c++14 -fopenmp Parallel_pca.cpp -o svd 
3. run parallel_pca as ./svd 784 784

- Steps to followed for serial code : 

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
