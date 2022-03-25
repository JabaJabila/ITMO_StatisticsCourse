clc;
clear;
pkg load statistics;


function stats = Kolmogorov_norm(a, sigma, n)
  i = [1 : n];
  x = sort(normrnd(a, sigma, 1, n));
  F = normcdf(x, a, sigma);
  stats = sqrt(n) * max(max(abs(F(i) - (i - 1) / n), (abs(F(i) - i / n))));
endfunction

function stats = Kolmogorov_unif(a, b, n)
  i = [1 : n];
  x = sort(unifrnd(a, b, 1, n));
  F = unifcdf(x, a, b);
  stats = sqrt(n) * max(max(abs(F(i) - (i - 1) / n), max(abs(F(i) - i / n))));
endfunction

function stats = Smirnov_norm(a, sigma, n)
  i = [1 : n];
  x = sort(normrnd(a, sigma, 1, n));
  F = normcdf(x, a, sigma);
  stats = 1/(12 * n) + sum((F(i) - (2 * i - 1) / (2 * n)).^ 2);
endfunction

function stats = Smirnov_unif(a, b, n)
  i = [1 : n];
  x = sort(unifrnd(a, b, 1, n));
  F = unifcdf(x, a, b);
  stats = 1/(12 * n) + sum((F(i) - (2 * i - 1) / (2 * n)).^ 2);
endfunction

function CheckStats(kolm, smirn)
    if (kolm < 1.224)
    printf("gamma = 0.9 Kolmogorov - OK!\n")
  else
    printf("gamma = 0.9 Kolmogorov - ERROR!\n")
  endif
  if (kolm < 1.358)
    printf("gamma = 0.95 Kolmogorov - OK!\n")
  else
    printf("gamma = 0.9 Kolmogorov - ERROR!\n")
  endif
  if (smirn < 0.35)
    printf("gamma = 0.9 Smirnov - OK!\n")
  else
    printf("gamma = 0.9 Smirnov - ERROR!\n")
  endif
  if (smirn < 0.46)
    printf("gamma = 0.95 Smirnov - OK!\n")
  else
    printf("gamma = 0.95 Smirnov - ERROR!\n")
  endif
endfunction

n1 = 10^4;
n2 = 10^6;

a = 0;
sigma = 4;
b = 3;

kolm = Kolmogorov_norm(a, sigma, n1);
smirn = Smirnov_norm(a, sigma, n1);

printf("\nNormal: (n=%d) Kolmogorov = %d; Smirnov = %d\n\n", 
n1, kolm, smirn);
CheckStats(kolm, smirn);

kolm = Kolmogorov_norm(a, sigma, n2);
smirn = Smirnov_norm(a, sigma, n2);

printf("\nNormal: (n=%d) Kolmogorov = %d; Smirnov = %d\n\n", 
n2, kolm, smirn);
CheckStats(kolm, smirn);

kolm = Kolmogorov_unif(a, b, n1);
smirn = Smirnov_unif(a, b, n1);

printf("\nUniform: (n=%d) Kolmogorov = %d; Smirnov = %d\n\n", 
n1, kolm, smirn);
CheckStats(kolm, smirn);

kolm = Kolmogorov_unif(a, b, n2);
smirn = Smirnov_unif(a, b, n2);

printf("\nUniform: (n=%d) Kolmogorov = %d; Smirnov = %d\n\n", 
n2, kolm, smirn);
CheckStats(kolm, smirn);
