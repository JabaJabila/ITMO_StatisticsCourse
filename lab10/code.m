clc;
clear;
pkg load statistics;

N = 10
p = 0.65
q = 1 - p

vec_p = p * ones(N - 2, 1);
vec_q = q * ones(N - 2, 1);
A = diag(vec_p, 1);
B = diag(vec_q, -1);
P = A + B;

P(1, 1) = 1; P(N, N) = 1;
P(1, 2) = 0; P(N-1, N) = p;

printf('Матрица:\n');
disp(P);

p01 = rand(1, N - 4);
p02 = p01 / sum(p01);
p0 = [0, 0, p02, 0, 0]

x = zeros(1, N);
x(1) = p0(1);
for i = 2 : N
  x(i) = x(i - 1) + p0(i);
end

u = rand(1, 1);
start = 1;

while (u > x(start))
  start = start + 1;
end


if (p != q)
  px0 = ((q / p) ^ (start - 1) - (q / p) ^ (N - 1)) / (1 - (q / p) ^ (N - 1));
else
  px0 = 1 - (start - 1) / (N - 1);
endif

pxN = 1 - px0;

printf('Теоретичесие предельные вероятности:\n')
[px0, pxN]

j = start;
x0 = zeros(1, N);
x0(j) = 1;
m = 200;
px_prac = x0 * (P ^ m);
printf('\nПрактические предельные вероятности:\n');
[px_prac(1), px_prac(N)]

i = 1;
for k = 1 : m + 1
  Pn = x0 * P ^ k;
  Pn1(i, :) = Pn;
  i++;
end

figure(1);
plot(0 : m, Pn1)
set(gca, "fontsize", 30);
title('График изменения вероятности')
grid;


figure(2);
s(1) = j;
K = 100;
for k = 2 : K
  u = rand(1, 1);
  if (s(k - 1) == 1 || s(k - 1) == N )
    s(k) = s(k - 1);
  elseif (u < q)
    s(k) = s(k - 1) - 1;
  elseif
    s(k) = s(k - 1) + 1;
  endif
end

plot(s, 'r--*');
set(gca, "fontsize", 30);
title('Моделирование случайного блуждания с поглощением')