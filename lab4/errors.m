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
  result = hi >= chi2inv(gamma, k - 3);
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
  result = hi >= chi2inv(gamma, k - 3);
endfunction

function norm_first_errors(mu, sigma, dmu, dsigma, n, gamma, count)
  errors = 0;
  for i = 1 : count
    errors += norm_criteria(mu, sigma, dmu, dsigma, n, gamma);
  endfor
  printf("%g\n", errors / count);
endfunction

function norm_second_errors(mu, sigma, dmu, dsigma, n, gamma, count)
  errors = 0;
  for i = 1 : count
    errors += norm_criteria(mu, sigma, dmu, dsigma, n, gamma);
  endfor
  printf("%g\n", (count - errors) / count);
endfunction

function unif_first_errors(a, b, da, db, n, gamma, count)
  errors = 0;
  for i = 1 : count
    errors += unif_criteria(a, b, da, db, n, gamma);
  endfor
  printf("%g\n", errors / count);
endfunction

function unif_second_errors(a, b, da, db, n, gamma, count)
  errors = 0;
  for i = 1 : count
    errors += unif_criteria(a, b, da, db, n, gamma);
  endfor
  printf("%g\n", (count - errors) / count);
endfunction

printf("Вер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^4, 0.9)
norm_first_errors(1, 3, 0, 0, 10^4, 0.9, 100);
##printf("Вер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^4, 0.95)
##norm_first_errors(1, 3, 0, 0, 10^4, 0.95, 100);
##printf("Вер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^4, 0.99)
##norm_first_errors(1, 3, 0, 0, 10^4, 0.99, 100);
##
##printf("\nВер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^6, 0.9)
##norm_first_errors(1, 3, 0, 0, 10^6, 0.9, 100);
##printf("Вер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^6, 0.95)
##norm_first_errors(1, 3, 0, 0, 10^6, 0.95, 100);
##printf("Вер. ошибки первого рода для нормального распределения для n = %d, gamma = %g:\t", 10^6, 0.99)
##norm_first_errors(1, 3, 0, 0, 10^6, 0.99, 100);
##
printf("\nВер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^4, 0.9)
unif_first_errors(2, 5, 0, 0, 10^4, 0.9, 100);
##printf("Вер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^4, 0.95)
##unif_first_errors(2, 5, 0, 0, 10^4, 0.95, 100);
##printf("Вер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^4, 0.99)
##unif_first_errors(2, 5, 0, 0, 10^4, 0.99, 100);
##
##printf("\nВер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^6, 0.9)
##unif_first_errors(2, 5, 0, 0, 10^6, 0.9, 100);
##printf("Вер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^6, 0.95)
##unif_first_errors(2, 5, 0, 0, 10^6, 0.95, 100);
##printf("Вер. ошибки первого рода для равномерного распределения для n = %d, gamma = %g:\t", 10^6, 0.99)
##unif_first_errors(2, 5, 0, 0, 10^6, 0.99, 100);


printf("\n\nВер. ошибки второго рода для нормального распределения для n = %d, gamma = %g (сдвиг по обоим параметрам = %g):\t", 10^6, 0.95, 0.008)
norm_second_errors(1, 3, 0.008, 0.008, 10^6, 0.95, 100);
##printf("Вер. ошибки второго рода для нормального распределения для n = %d, gamma = %g (сдвиг по обоим параметрам = %g):\t", 10^6, 0.95, 0.012)
##norm_second_errors(1, 3, 0.012, 0.012, 10^6, 0.95, 100);
##printf("Вер. ошибки второго рода для нормального распределения для n = %d, gamma = %g (сдвиг по обоим параметрам = %g):\t", 10^6, 0.95, 0.015)
##norm_second_errors(1, 3, 0.015, 0.015, 10^6, 0.95, 100);

printf("\nВер. ошибки второго рода для равномерного распределения для n = %d, gamma = %g (сдвиг по первому параметру = %g):\t", 10^6, 0.95, 0.001)
unif_second_errors(2, 5, 0.001, 0, 10^6, 0.95, 100);
##printf("Вер. ошибки второго рода для равномерного распределения для n = %d, gamma = %g (сдвиг по первому параметру = %g):\t", 10^6, 0.95, 0.002)
##unif_second_errors(2, 5, 0.002, 0, 10^6, 0.95, 100);
##printf("Вер. ошибки второго рода для равномерного распределения для n = %d, gamma = %g (сдвиг по первому параметру = %g):\t", 10^6, 0.95, 0.004)
##unif_second_errors(2, 5, 0.004, 0, 10^6, 0.95, 100);
