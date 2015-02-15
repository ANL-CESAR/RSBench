echo "" > sweepOutput

for innerDim in 16 32 64 128 256 512 1024; do
    for batchDim in 1 2 4 8 16 32 64 128 256 512 1024; do
        echo "Inner Dim [$innerDim], Batch Dim [$batchDim]" >> sweepOutput
        ./rsbench -s small -d -i $innerDim -b $batchDim | grep Runtime >> sweepOutput
    done
done