function [x0,x1,k]=seccz(x0,x1,f,erro)
  
  %EXEMPLO: [x0,x1,k]=seccz(0,1,f,0.0005)
  % f= @(x) x .^ 3 - 9 * x + 3
  
  disp('-----------------')
  disp('Método da Secante')
  disp('-----------------')
  
  disp("  k          x0         f(x0)        x1        f(x1)     abs(x0-p)")
  
  %disp('Insira o valor para x0:') %0
  %x0 = input('')
  
  %disp('Insira o valor para x1:') %1
  %x1 = input('')
  
  %disp('Defina a função usando o @(x):')
  %disp('EXEMPLO: @(x) x .^ 3 - 9 * x + 3')
  %f = input('')
  
  %disp('Defina o erro para o critério de parada:') %0.0005
  %erro = input('')
  
  k=0;
  while k<100
    p=(x0*f(x1)-x1*f(x0))/(f(x1)-f(x0));
    if f(p)==0 || abs(x0-p)<erro
      break
    else 
      
      printf('%f   %f   %f   %f   %f   %f\n', k,x0,f(x0),x1,f(x1),abs(x0-p))
      
      k=k+1;
      x0=x1;
      x1=p;
    endif
  endwhile
  disp('-----------------')
endfunction