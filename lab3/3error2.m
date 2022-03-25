clc;
clear;
pkg load statistics;


function stats = Kolmogorov_norm(a, da, sigma, dsigma, n)
  i = [1 : n];
  x = sort(normrnd(a, sigma, 1, n));
  F = normcdf(x, a + da, sigma + dsigma);
  stats = sqrt(n) * max(max(abs(F(i) - (i - 1) / n), (abs(F(i) - i / n))));
endfunction

function stats = Kolmogorov_unif(a, da, b, db, n)
  i = [1 : n];
  x = sort(unifrnd(a, b, 1, n));
  F = unifcdf(x, a + da, b + db);
  stats = sqrt(n) * max(max(abs(F(i) - (i - 1) / n), max(abs(F(i) - i / n))));
endfunction

function stats = Smirnov_norm(a, da, sigma, dsigma, n)
  i = [1 : n];
  x = sort(normrnd(a, sigma, 1, n));
  F = normcdf(x, a + da, sigma + dsigma);
  stats = 1/(12 * n) + sum((F(i) - (2 * i - 1) / (2 * n)).^ 2);
endfunction

function stats = Smirnov_unif(a, da, b, db, n)
  i = [1 : n];
  x = sort(unifrnd(a, b, 1, n));
  F = unifcdf(x, a + da, b + db);
  stats = 1/(12 * n) + sum((F(i) - (2 * i - 1) / (2 * n)).^ 2);
endfunction

n = 10^4;

a = 0;
sigma = 4;
b = 3;

function error_2(a, da, b, db, sigma, dsigma, n, m)
  
  norm_dn1 = 0;
  norm_wn1 = 0;
  unif_dn1 = 0;
  unif_wn1 = 0;
  norm_dn2 = 0;
  norm_wn2 = 0;
  unif_dn2 = 0;
  unif_wn2 = 0;
  
  gamma1 = 0.9;
  dn1 = 1.224;
  wn1 = 0.35;
  
  gamma2 = 0.95;
  dn2 = 1.358;
  wn2 = 0.46;
  
  
  for (i = [1 : m])
    
    k_norm = Kolmogorov_norm(a, da, sigma, dsigma, n);
    w_norm = Smirnov_norm(a, da, sigma, dsigma, n);
    
    if (k_norm < dn1)
      norm_dn1 += 1;
    endif
    
    if (w_norm < wn1)
      norm_wn1 += 1;
    endif
    
    if (k_norm < dn2)
      norm_dn2 += 1;
    endif
    
    if (w_norm < wn2)
      norm_wn2 += 1;
    endif
    
    k_unif = Kolmogorov_unif(a, da, b, db, n);
    w_unif = Smirnov_unif(a, da, b, db, n);
    
    if (k_unif < dn1)
      unif_dn1 += 1;
    endif
    
    if (w_unif < wn1)
      unif_wn1 += 1;
    endif
    
    if (k_unif < dn2)
      unif_dn2 += 1;
    endif
    
    if (w_unif < wn2)
      unif_wn2 += 1;
    endif
    
  endfor
  
  printf("\ngamma = %d (alpha = %d); n = %d; m = %d\n", gamma1, 1 - gamma1, n, m);
  printf("Normal(Fn(a=%d,sigma=%d); Fx(a=%d,sigma=%d)): Kolmogorov error = %d \tSmirnov error = %d\n", a, sigma, a + da, sigma + dsigma, norm_dn1 / m, norm_wn1 / m);
  
  printf("\ngamma = %d (alpha = %d); n = %d; m = %d\n", gamma2, 1 - gamma2, n, m);
  printf("Normal(Fn(a=%d,sigma=%d); Fx(a=%d,sigma=%d)): Kolmogorov error = %d \tSmirnov error = %d\n", a, sigma, a + da, sigma + dsigma, norm_dn2 / m, norm_wn2 / m);
  
  printf("\ngamma = %d (alpha = %d); n = %d; m = %d\n", gamma1, 1 - gamma1, n, m);
  printf("Uniform(Fn(a=%d,b=%d); Fx(a=%d,b=%d)): Kolmogorov error = %d \tSmirnov error = %d\n", a, b, a + da, b + db, unif_dn1 / m, unif_wn1 / m);
  
  printf("\ngamma = %d (alpha = %d); n = %d; m = %d\n", gamma2, 1 - gamma2, n, m);
  printf("Normal(Fn(a=%d,b=%d); Fx(a=%d,b=%d)): Kolmogorov error = %d \tSmirnov error = %d\n", a, b, a + da, b + db, unif_dn2 / m, unif_wn2 / m);
endfunction

error_2(a, 0.005, sigma, 0.005, b, 0.005, n, 100);
error_2(a, 0.01, sigma, 0.01, b, 0.01, n, 100);
error_2(a, 0.02, sigma, 0.02, b, 0.02, n, 100);
error_2(a, 0.03, sigma, 0.03, b, 0.03, n, 100);
error_2(a, 0.03, sigma, 0.03, b, 0.03, n, 100);
error_2(a, 0.05, sigma, 0.05, b, 0.05, n, 100);
