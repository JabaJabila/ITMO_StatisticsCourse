clc
clear
pkg load statistics

a = 1;
b = 8;
gamma = 0.95;

n1 = 10^4;

T = norminv((gamma + 1) / 2);
x = unifrnd(a, b, n1, 1);
Y = (b - a) * (1 ./ (sqrt(1 + x .^ 3)));
I = mean(Y)
d = T * std(Y) / sqrt(n1)
gov_int1 = [I - d, I + d]

n2 = 10^6;

x = unifrnd(a, b, n2, 1);
Y = (b - a) * (1 ./ (sqrt(1 + x .^ 3)));
I = mean(Y)
d = T * std(Y) / sqrt(n2)
dov_int2 = [I - d, I + d]

Ireal = quad('1 ./ (sqrt(1 + x .^ 3))', 1, 8)