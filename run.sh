#!/bin/bash
rm output.txt
bash compile.sh
echo -ne "SERIAL"
echo "...................." >> output.txt
echo "...................." >> output.txt
./serial_pca 784 784 8.txt 90 >> output.txt
echo -ne "PARALLEL"
echo "...................." >> output.txt
echo "...................." >> output.txt
echo "PARALLEL PART" >> output.txt
for i in 1 2 4 6 8 16;
do
	echo -ne "."
	echo "...................." >> output.txt
    echo "...................." >> output.txt
	echo "Results for thread count = ${i}" >> output.txt
	./parallel_pca 784 784 ${i} 8.txt 90 >> output.txt
	echo "...................." >> output.txt
    echo "...................." >> output.txt
	echo -ne "."
done

echo "Completed"
