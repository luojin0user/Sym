tic
[row_idx, col_idx] = ndgrid(1:600, 1:600);
toc;
D = gpuDevice;
wait(D)
tic
row_idx = gpuArray(row_idx);
col_idx = gpuArray(col_idx);
wait(D)
toc;

