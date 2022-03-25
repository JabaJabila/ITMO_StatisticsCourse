clc;
clear;
pkg load statistics;

m = 10 ^ 2;
n = 10 ^ 4;
gamma = 0.95;
a = 1;
b = 7;
t0 = 3.2;

T = norminv((1 + gamma) / 2);
Fx0 = unifcdf(t0, a, b)

lower = [];
upper = [];

for i = 1 : m
  X = unifrnd(a, b, n, 1);
  Fx = mean(X < t0);
  d = T * sqrt(Fx * (1 - Fx)) / sqrt(n);
  lower = [lower, Fx - d];
  upper = [upper, Fx + d];
endfor

figure(1)
plot(1 : m, lower, 1 : m, upper, 1 : m, linspace(Fx0, Fx0, m));
title("U(1, 7) (gamma = 0.95)", "fontsize", 24);

k = 100;

for gamma = [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
  count = [];
  T = norminv((1 + gamma) / 2);
  
  for i = 1 : k
    lower = [];
    upper = [];
    
    for j = 1 : m
      X = unifrnd(a, b, n, 1);
      Fx = mean(X < t0);
      d = T * sqrt(Fx * (1 - Fx)) / sqrt(n);
      lower = [lower, Fx - d];
      upper = [upper, Fx + d];
    endfor
    
    count = [count, sum(lower > Fx0) + sum(upper < Fx0)];
  endfor
  
  printf("gamma = %g: %g (1 - gamma = %g)\n", gamma, mean(count), 1 - gamma);
endfor

printf("----------")

