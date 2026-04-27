function [X, R] = imstack2vectors(S, MASK)
[M, N, n] = size(S);
if nargin == 1
   MASK = true(M, N);    %%  M*N 个1的逻辑矩阵
else
   MASK = MASK ~= 0;
end

[I, J] = find(MASK);    %%返回非0元素
R = [I, J];             %%非0元素的坐标矩阵

Q = M*N;
X = reshape(S, Q, n);

MASK = reshape(MASK, Q, 1);

X = X(MASK, :);

