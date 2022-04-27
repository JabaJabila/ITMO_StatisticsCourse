clc;
clear;
pkg load statistics;


mu = 1;
sigma = 3;
n = 10^6;


x = [-15 : 0.01 : 15];
y = normpdf(x, mu, sigma);

d = sort(normrnd(mu, sigma, n, 1));
histogram = hist(d, 120);
h = (d(end) - d(1)) / 120;
xx = d(1) : h : d(end) - h;
yy = histogram / n / h;
[a, b] = stairs(xx, yy);

figure(1)
plot(x, y, a, b), grid
set(gca, "fontsize", 24);
title(strcat("N (", num2str(mu),",",num2str(sigma),")"));

a1 = 2
b1 = 5

x = [a1 - 0.2 : 0.01 : b1 + 0.2];
y = unifpdf(x, a1, b1);
d = sort(unifrnd(a1, b1, 1, n));
histogram = hist(d, 120);
h = (d(end) - d(1)) / 120;
xx = d(1) : h : d(end) - h;
yy = histogram / n / h;
[a, b] = stairs(xx, yy);

figure(2)
plot(x, y, a, b), grid
set(gca, "fontsize", 24);
title(strcat("U (", num2str(a1),",",num2str(b1),")"));