clc;
clear;
pkg load statistics;

function result = norm_criteria(mu, sigma, dmu, dsigma, n, gamma)
  X = sort(normrnd(mu, sigma, 1, n));
  q1 = sum(X) / n;
  q2 = sqrt(sum((X - mu).^2) / n);
 
  k = ceil(n ^ (1 / 3));
  h = (X(n) - X(1)) / k;
  histogram = hist(X, k);

  q1 = q1 + dmu;
  q2 = q2 + dsigma;
  hi = 0;
  i = 1;
  m = 0;
  while (i <= k)
    oi = i - 1;
    tmp = histogram(i);
    while (tmp < 5 && i < k)
      i = i + 1;
      tmp = tmp + histogram(i);
    endwhile
    i_p = normcdf(X(1) + h * i, q1, q2) - normcdf(X(1) + h * oi, q1, q2);
    hi = hi + (tmp - n * i_p) ^ 2 / (n * i_p);
    i = i + 1;
    m = m + 1;
  endwhile
  
  printf("\nгамма = %g, порог = %g, степеней свободы = %g, статистика = %g\n",
          gamma, chi2inv(gamma, k - 3), k - 3, hi);
          
  result = "да\n";
  if (hi >= chi2inv(gamma, k - 3))
    result = "нет\n";
  endif
endfunction

function result = unif_criteria(a, b, da, db, n, gamma)
  X = sort(unifrnd(a, b, 1, n));
  q1 = X(1);
  q2 = X(n);
  
  k = ceil(n ^ (1 / 3));
  h = (X(n) - X(1)) / k;
  histogram = hist(X, k);

  q1 = q1 + da;
  q2 = q2 + db;
  hi = 0;
  i = 1;
  m = 0;
  while (i <= k)
    tmp = histogram(i);
    i_p = unifcdf(q1 + h * i, q1, q2) - unifcdf(q1 + h * (i - 1), q1, q2);
    hi = hi + (tmp - n * i_p)^2 / (n * i_p);
    i = i + 1;
    m = m + 1;
  endwhile
  
  printf("\nгамма = %g, порог = %g, степеней свободы = %g, статистика = %g\n",
          gamma, chi2inv(gamma, k - 3), k - 3, hi);
          
  result = "да\n";
  if (hi >= chi2inv(gamma, k - 3))
    result = "нет\n";
  endif
endfunction

mu = 1;
sigma = 3;
n = 10^6;

printf(strcat("\t\tN (", num2str(mu),",",num2str(sigma),")\n"));
res = norm_criteria(mu, sigma, 0, 0, n, 0.9);
printf("Принимаем: %s", res);
res = norm_criteria(mu, sigma, 0, 0, n, 0.95);
printf("Принимаем: %s", res);
res = norm_criteria(mu, sigma, 0, 0, n, 0.99);
printf("Принимаем: %s", res);

a = 2;
b = 5;

printf(strcat("\n\n\n\t\tU (", num2str(a),",",num2str(b),")\n"));
res = unif_criteria(a, b, 0, 0, n, 0.9);
printf("Принимаем: %s", res);
res = unif_criteria(a, b, 0, 0, n, 0.95);
printf("Принимаем: %s", res);
res = unif_criteria(a, b, 0, 0, n, 0.99);
printf("Принимаем: %s", res);