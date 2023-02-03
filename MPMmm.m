  disp('-------------------------------------------------------------------')
  disp('-------------------------------------------------------------------')
  disp('Insira a EDO:')
  disp('')
  edo = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o valor inicial de X, chamado de x0:')
  disp('')
  X0 = input('');
  disp('Insira o valor inicial de Y, chamado de y0:')
  disp('')
  y0 = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o passo que será utilizado, chamado por h:')
  disp('')
  h = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira até qual x você deseja chegar, chamado de xn:')
  disp('')
  Xn = input('');
  disp('-------------------------------------------------------------------')

  %%%FIM DO INPUT%%%
  disp("  i        x(i)       y(i)       k1(i)       y1/2(i)       k2(i)")

  n=(xn-x0)/h;

  x(1)=x0;
  y(1)=y0;

  for i = 1:(n+1)

    x(i+1) = x(i)+h;
    K1 = edo(x(i),y(i));

    ym = y(i)+K1*(h/2); %retira para euler%
    K2 = edo((x(i)+(h/2)),ym); %retira para euler%

    %y(i+1) = y(i)+K1*h %troca com a de baixo para euler%
    y(i+1) = y(i)+K2*h; %retira para euler%
    if i>n
      printf('%f   %f%   f\n', i, x(i),y(i))
    else
     printf('%f   %f   %f   %f   %f   %f\n', i, x(i),y(i), K1, ym, K2)
    endif
  endfor
  %%SAÍDA%%
