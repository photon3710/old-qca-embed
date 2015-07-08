#!/bin/bash

echo "changing to src directory"
cd ..

# assert input argument
if [ $# -eq 0 ]; then
	echo "Need to supply the loop count"
	exit
fi

echo "Running gen_test the given number of times"
for (( i=0;i<$1;i++)); do
	python gen_test.py $2
done

echo "changing back to sh-scipts directory"
cd ./sh-scripts
