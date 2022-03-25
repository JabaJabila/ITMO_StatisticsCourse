clc
clear
pkg load statistics

a = 2
k = 3 
c = 4.3 
gamma = 0.95 

n1 = 10^4
T = norminv((1 + gamma) / 2);

x = rand(n1, k); 
f = a .^ x; 
z = sum(f'); 
v1 = sum(z <= c) / n1
diff1 = T * sqrt(v1 * (1 - v1) / n1)
dov_int1 = [v1 - diff1, v1 + diff1]

n2 = 10^6
x = rand(n2, k); 
f = a.^x;
z = sum(f'); 
v2 = sum(z <= c) / n2
diff2 = T * sqrt(v2 * (1 - v2) / n2) 
dov_int2 = [v2 - diff2, v2 + diff2]