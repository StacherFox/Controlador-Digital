clc;

kmax = 20;

k = 0:kmax;

r = 0:Ts:kmax*Ts;
e = zeros(1,length(k));
c = zeros(1,length(k));
y = zeros(1,length(k));

n = 1;
y(n) = 0;
e(n) = r(n) - y(n);
c(n) = numc(1)*e(n);

n = 2;
y(n) = numgz(2)*c(n-1) - dengz(2)*y(n-1);
e(n) = r(n) - y(n);
c(n) = numc(1)*e(n) + numc(2)*e(n-1) - denc(2)*c(n-1);

n = 3;
y(n) = numgz(2)*c(n-1) + numgz(3)*c(n-2) - dengz(2)*y(n-1) - dengz(3)*y(n-2);
e(n) = r(n) - y(n);
c(n) = numc(1)*e(n) + numc(2)*e(n-1) + numc(3)*e(n-2) - denc(2)*c(n-1) - denc(3)*c(n-2);

for n = 4:length(k)
   y(n) = numgz(2)*c(n-1) + numgz(3)*c(n-2) + numgz(4)*c(n-3) - dengz(2)*y(n-1) - dengz(3)*y(n-2) - dengz(4)*y(n-3);
   e(n) = r(n) - y(n);
   c(n) = numc(1)*e(n) + numc(2)*e(n-1) + numc(3)*e(n-2) - denc(2)*c(n-1) - denc(3)*c(n-2);
end

plot(k*Ts, y, '*');