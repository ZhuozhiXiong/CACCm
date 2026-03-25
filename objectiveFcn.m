function f = objectiveFcn(optimInput)
x1 = optimInput(1);
x2 = optimInput(2);
x3 = optimInput(3);
x4 = 0;
k1 = optimInput(4);
k2 = optimInput(5);
h = 1.2;
theta = 0;
m = 3;
n = 1;

lb=[-0.031 -0.062 0.75 0.25];
ub=[0.431 0.862 2.25 1.75];
N = 5000;

HDV_v1 = [x1 x2 x3 x4 k1 k2 h theta m n];
HDV_v2 = repmat(HDV_v1,N,1);
HDV_v = num2cell(HDV_v2, 1);
V = prod(abs(ub-lb));
n = length(lb);
rng(1)
sample = (ub-lb) .* rand(N, n) + repmat(lb,N,1);
sample_arguments = num2cell(sample, 1);
results = cell2mat(arrayfun(@ obj_function, sample_arguments{:},HDV_v{:}, 'UniformOutput', 0));
f = sum(results) * V / N;
end