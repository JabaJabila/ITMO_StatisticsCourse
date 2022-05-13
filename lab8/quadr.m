clc
clear
pkg load statistics

x_min = 0;
x_max = 2;
n = 80;
s = 2.3;

ai = [2.1, -3.2, 1.4];

printf("Функция: y = %gx^2 + %gx + %g\tn = %d\t s = %g\n\n", ai(3), ai(2), ai(1), n, s);

X = linspace(x_min, x_max, n)';
A = [];
i = [1, 2, 3];

A(:, i) = X .^ (i - 1);
y = A * ai';
Z = s * randn(n, 1);
Y = y + Z;

ai_function = polyfit(X, Y, 2);
ai_matrix = (A' * A) ^ -1 * A' * Y;


printf("Параметры:  \t\t\t %f \t %f \t %f\n", ai(3), ai(2), ai(1));
printf("Расчёт встроенной функцией:  \t %f \t %f \t %f\n",  ai_function(1), ai_function(2), ai_function(3));
printf("Матричный расчёт:  \t\t %f \t %f \t %f\n\n", ai_matrix(3), ai_matrix(2), ai_matrix(1));

Y_matrix = A * ai_matrix;
Y_function = polyval(ai_function, X);

r = Y_matrix - Y;
printf("Ортогональность: %d\n", r' * Y_matrix);
s_n = sqrt(r' * r / (n - 3));
printf("Расчёт уровня шумов = %d\n", s_n);

plot(X, Y, '.', X, Y_function, '+', X, Y_matrix, 'o', X, y, '-');
legend("Точки", "Встроенная функция", "Матричный расчёт", "Функция");
axis("tight");

