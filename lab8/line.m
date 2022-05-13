clc
clear
pkg load statistics

x_min = 0;
x_max = 2;
n = 80;
s = 2.3;

ci = [-2.1, 1.8];

printf("Функция: y = %gx + %g\tn = %d\t s = %g\n\n", ci(2), ci(1), n, s);

X = linspace(x_min, x_max, n)';
A = [];
i = [1, 2];

A(:, i) = X .^ (i - 1);
y = A * ci';
Z = s * randn(n, 1);
Y = y + Z;

ci_function = polyfit(X, Y, 1);
ci_matrix = (A' * A) ^ -1 * A' * Y;

cov_xy = X' * Y / n - mean(X) * mean(Y);
var_x = X' * X / n - (mean(X)) ^ 2;
cov_matrix = [mean(Y) - cov_xy / var_x * mean(X); cov_xy / var_x];

printf("Параметры:  \t\t\t %f \t %f\n", ci(2), ci(1));
printf("Расчёт встроенной функцией:  \t %f \t %f\n",  ci_function(1), ci_function(2));
printf("Матричный расчёт:  \t\t %f \t %f\n", ci_matrix(2), ci_matrix(1));
printf("Ковариация:  \t\t\t %f \t %f\n\n", cov_matrix(2), cov_matrix(1));

Y_matrix = A * ci_matrix;
Y_function = polyval(ci_function, X);
Y_cov_matrix = A * cov_matrix;

r = Y_matrix - Y;
printf("Ортогональность: %d\n", r' * Y_matrix);
s_n = sqrt(r' * r / (n - 3));
printf("Расчёт уровня шумов = %d\n", s_n);


plot(X, Y, '.', X, Y_function, '+', X, Y_matrix, 'o', X, y, '-', X, Y_cov_matrix, '^');
legend("Точки", "Встроенная функция", "Матричный расчёт", "Функция", "Ковариация");
axis("tight");