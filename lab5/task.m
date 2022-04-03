clc;
clear;
pkg load statistics;

function printTable(name, arg1, arg2, theory, experiment, n, m)
  theory = sqrt(theory);
  printf("\t\t\t\t%s(%d, %d)\n\n", name, arg1, arg2)
  printf("n = %d / m = %d\t|", n, m)
  printf("\tmean(x)\t\t|\tmed(x)\t\t|\t(x(1)+x(n))/2\n")
  printf("------------------------|-----------------------|-----------------------|-----------------------\n")
  printf("sigma-теор.: \t\t|\t%f\t|\t%f\t|\t%f\n", theory(1), theory(2), theory(3))
  printf("------------------------|-----------------------|-----------------------|-----------------------\n")
  printf("sigma-прак.: \t\t|\t%f\t|\t%f\t|\t%f\n", experiment(1), experiment(2), experiment(3))
  printf("------------------------|-----------------------|-----------------------|-----------------------\n")
  printf("Разница: \t\t|\t%f\t|\t%f\t|\t%f\n\n\n", abs(theory(1) - experiment(1)), abs(theory(2) - experiment(2)), abs(theory(3) - experiment(3)))
endfunction

function answer = theoryNormal(sigma, n)
  answer = [(sigma ^ 2) / n, (pi * sigma ^ 2) / (2 * n), (0.4 * sigma ^ 2) / log(n)];
endfunction

function answer = theoryLaplas(u, n)
  answer = [(2 * u ^ 2) / n, (u ^ 2) / n, (0.9 * u ^ 2)];
endfunction

function answer = theoryUniform(delta, n)
  answer = [(delta ^ 2) / (12 * n), (delta ^ 2) / (4 * n), (delta ^ 2) / (2 * n ^ 2)];
endfunction

function answer = experiment(X)
  answer = [std(mean(X)), std(median(X)), std((max(X) + min(X)) / 2)];
endfunction

function do3Distributions(a1, sigma, a2, delta, a3, u, n, m)
  X = normrnd(a1, sigma, n, m);
  printTable("Нормальное распределение N", a1, sigma, theoryNormal(sigma, n), experiment(X), n, m);

  X = unifrnd(a2 - delta / 2, a2 + delta / 2, n, m);
  printTable("Равномерное распределение U", a2 - delta / 2, a2 + delta / 2, theoryUniform(delta, n), experiment(X), n, m);

  X = a3 + exprnd(u, n, m) - exprnd(u, n, m);
  printTable("Распределение Лапласа L", a3, u, theoryLaplas(u, n), experiment(X), n, m);
endfunction

do3Distributions(2, 3, 3, 5, 3, 4, 100, 100)
do3Distributions(2, 3, 3, 5, 3, 4, 10000, 100)