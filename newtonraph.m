function[p,k]=newtonraph(x,f,erro,df)
       %%EXEMPLO%%
% x = 1.5
% f = @(x) x .^ 3 - 9 * x + 3
% erro = 0.003
% df = @(x) 3 * x .^ 2 + 9

% [p,k]=newtonraph(1.5,f,0.003,df)

disp('-----------------')
disp('Método de Newton-Raphson')
disp('-----------------')

disp("  k          x         f(x)          df(x)        p       abs((p-x)/p)")

k = 0;
while k<100
  
  p = x - (f(x)/df(x));
  
  if abs(f(p))<erro
    break
    % abs(p - x0)<erro
    % abs(p - x0/p)<erro
  endif
  
  printf('%f   %f   %f   %f   %f   %f\n', k,x,f(x),df(x),p,abs((p-x)/p))
  x = p;
  k = k+1;
endwhile
disp('-----------------')
endfunction