disp('-------------------------------------------------------------------')
disp('                   EDO pelo método de Ralson                       ')
disp('-------------------------------------------------------------------')

disp('Insira a EDO:')
disp('Exercício: @(x,y) 0.2*10*((100*(20-25)*(20-x))/100-2.5*x)+0*y')
disp('')
edo1 = input('');

disp('-------------------------------------------------------------------')
disp('Insira o valor inicial de X, chamado de x0:')
disp('Exercício: 0')
disp('')
x0 = input('');

disp('Insira o valor inicial de Y, chamado de y0:')
disp('Exercício: 0.7')
disp('')
y0 = input('');

%disp('Insira o valor inicial de Z, chamado de z0:')
%disp('Exercício: NÃO TINHA')
%disp('')
%z0 = input('');

disp('-------------------------------------------------------------------')
disp('Insira o passo que será utilizado, chamado por h:')
disp('Exercício: 0.1')
disp('')
h = input('');

disp('-------------------------------------------------------------------')
disp('Insira até qual x você deseja chegar, chamado de xn:')
disp('Exercício: 1')
disp('')
xn = input('');
disp('-------------------------------------------------------------------')
disp('                            RESULTADO                              ')
disp('-------------------------------------------------------------------')
%%%FIM DO INPUT%%%

  x(1)=x0;
  y(1)=y0;
  %z(1)=z0;
  
  n=(xn-x0)/h;
  
  %edo1=@(x,y,z) -0.5*y;
  %edo2=@(x,y,z) 4-0.3*z-0.1*y;
  
  disp("i\t    x(i)\t    y(i)\t   x(i+1)\t    y(i+1)\t");
  %disp("i\t x(i)\t y(i)\t z(i)\t k1y\t k1z\t k2y\t k2z\t x(i+1)\t y(i+1)\t z(i+1)"); FORMATO COM Z
  
  for i=1:n
    
    k1y=edo1(x(i),y(i));
    %k1y=edo1(x(i),y(i),z(i));  FORMATO COM Z
    %k1z=edo2(x(i),y(i),z(i));
    
    xapoi=x(i)+(3/4)*h;
    yapoi=y(i)+(3/4)*k1y*h;
    %zapoi=z(i)+(3/4)*k1z*h;
    
    k2y=edo1(xapoi,yapoi);
    %k2y=edo1(xapoi,yapoi,zapoi); FORMATO COM Z 
    %k2z=edo2(xapoi,yapoi,zapoi);
    
    x(i+1)=x(i)+h;
    y(i+1)=y(i)+((1/3)*k1y+(2/3)*k2y)*h; 
    %z(i+1)=z(i)+((1/3)*k1z+(2/3)*k2z)*h; 
 
    printf("%i\t %f\t %f\t %f\t %f\t \n",i,x(i),y(i),x(i+1),y(i+1));
    %printf("%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",i,x(i),y(i),z(i),x(i+1),y(i+1),z(i+1)); FORMATO COM Z
    
  endfor

  printf("%i\t %f\t %f\t \n",i+1,x(i+1),y(i+1));
  %printf("%i\t%.4f\t%.4f\t%.4f\n",i+1,x(i+1),y(i+1),z(i+1)); FORMATO COM Z
  
  disp('-------------------------------------------------------------------')