function[x,y] = euler3ordem(edo,x0,y0,h,xn)
n=(xn - x0)/h;
x(1) = x0;
y(1) = y0;
for i = 1:n
  x(i+1) = x(i) + h;
  K1 = edo(x(i),y(i));
  K2 = edo(x(i) + (h/2), y(i) + K1 * (h/2));
  K3 = edo(x(i) + h, y(i) - (K1 * h) + (2 * K2 * h));
  y(i+1) = y(i) + (K1 + 4*K2 + K3) * (h/6)
endfor
endfunction
