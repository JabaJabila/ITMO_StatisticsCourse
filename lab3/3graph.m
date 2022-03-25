clc;
clear;
pkg load statistics;


a = 0;
sigma = 4;
n = 100;
gamma = 0.95;
u_gamma = 1.36;

x_p = [-10 : 1/n : 10];
y_p = normcdf(x_p, a, sigma);

x = sort(normrnd(a, sigma, n, 1));
y = [1/n: 1/n : 1];
[x_ps, y_ps] = stairs(x, y);

lower = max(0, y_ps - u_gamma / sqrt(n));
upper = min(1, y_ps + u_gamma / sqrt(n));

figure(1);
plot(x_p, y_p, x_ps, y_ps, x_ps, lower, x_ps, upper);
set(gca, "linewidth", 1, "fontsize", 18);
title("Normal(a = 0, sigma = 4) (gamma = 0.95)", "fontsize", 24);


a = 0;
b = 3;

x_p = [a : 1/n : b];
y_p = unifcdf(x_p, a, b);

x = sort(unifrnd(a, b, n, 1));
y = [1/n: 1/n : 1];
[x_ps, y_ps] = stairs(x, y);

lower = max(0, y_ps - u_gamma / sqrt(n));
upper = min(1, y_ps + u_gamma / sqrt(n));

figure(2);
plot(x_p, y_p, x_ps, y_ps, x_ps, lower, x_ps, upper);
set(gca, "linewidth", 1, "fontsize", 18);
title("Uniform(a = 0, b = 3) (gamma = 0.95)", "fontsize", 24);
