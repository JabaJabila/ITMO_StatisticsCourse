clc
clear
pkg load statistics

a = -1;
sigma = 1 / sqrt(2);
gamma = 0.95;

n1 = 10^4;

T = norminv((1 + gamma) / 2);
x = normrnd(a, sigma, n1, 1);
Y = sigma * sqrt(2 * pi) * sin(x);
I = mean(Y)
dI = T * std(Y) / sqrt(n1)
dov_int1 = [I - dI, I + dI]

n2 = 10^6;
x = normrnd(a, sigma, n2, 1);
Y = sigma * sqrt(2 * pi) * sin(x);
I = mean(Y)
dI = T * std(Y) / sqrt(n2)
dov_int2 = [I - dI, I + dI]

Ireal = quad('sin(x).*exp(-(x + 1) .^ 2)', -inf, inf)