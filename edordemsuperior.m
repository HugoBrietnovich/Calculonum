function[x,y,z] = euler2edos(edoy,edoz,x0,y0,z0,h,xn)
n=(xn - x0)/h;
x(1) = x0;
y(1) = y0;
z(1) = z0;
for i = 1:n
  x(i+1) = x(i) + h;
  K1y = edoy(x(i),y(i),z(i));
  K1z = edoz(x(i),y(i),z(i));
  y(i+1) = y(i) + K1y * h;
  z(i+1) = z(i) + K1z * h;
endfor
endfunction
