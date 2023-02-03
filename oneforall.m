%ONEFORALL%
%Versão 7.0Kb%
%CRÉDITOS AO CRIADOR%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Bem vindo ao One for All!!')
disp('  Este é um programa feito no Octave para rodar e servir de uso na ')
disp('matéria de Cálculo Numérico ao compilar todas as implementações dos')
disp('métodos numéricos utilizados.')
disp('')
disp('CRÉDITOS AO CRIADOR: Victor Hugo @hyuug4rts (AH EH!!)')

%%Para desabilitar esse pedido de formato, coloque "%" na frente de todas as linhas abaixo

disp('-------------------------------------------------------------------')
disp('Primeiro defina qual formato deseja utilizar:')
disp('')
disp('1) Format Long')
disp('Neste formato seu octave irá mostrar uma MAIOR quantidade de casas decimais.')
disp('')
disp('2) Format Short')
disp('Neste formato seu octave irá mostrar uma MENOR quantidade de casas decimais.')
forma = input('');

if forma == 1
  format long
  disp('Formato de representação LONGA em uso...')
endif
if forma == 2
  format short
  disp('Formato de representação CURTA em uso...')
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% até aqui

disp('-------------------------------------------------------------------')
disp('Escolha o que deseja fazer:')
opcao = menu('Escolha o método:', '1)MÉTODO: Bissecção', '2)MÉTODO: Falsa Posição', '3)MÉTODO: Newton-Raphson', '4)MÉTODO: Secante', '5)MÉTODO Gauss-Seidel', '6)MÉTODO Gauss-Jacobi', '7)Fatoração LU', '8)Métodos de EDO', '9)Interpolação Polinomial', '10)Integração Numérica','11)Ajuste de Curvas', '12)Plotar Gráficos', '13)Informações adicionais')

switch(opcao)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1 %MÉTODO Bissecção

%%%INPUT%%%

disp('-------------------------------------------------------------------')
disp('                     MÉTODO EM USO: Bissecção                      ')
disp('-------------------------------------------------------------------')
disp('Defina os intervalos do "chute inicial"')

disp('Insira o valor de a:')
disp('EX: 3.111')
disp('')
a = input('');

disp('Insira o valor de b:')
disp('EX: 3.222')
disp('')
b = input('');
disp('-------------------------------------------------------------------')
disp('Me informe a função que deseja aplicar o método da bissecção:')
disp('EX: @(x) sin(x)')
disp('')
f = input('');
disp('-------------------------------------------------------------------')
disp('Defina o erro a ser usado no critério de parada:')
disp('Ex: 0.001')
disp('')
erro = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
disp("                        TABELA DE ITERAÇÕES                        ")
disp('-------------------------------------------------------------------')
disp("  k        a         b          f(a)        f(b)       p       f(p)")

if f(a)*f(b) > 0
   error('Ta errado irmão, troca issae sobrinho ~By Bill!!')
else  
  k=0;
  
  while k<100 
        p = (a+b)/2;
 
    if abs(b-a)<erro %Critério de parada%
      format long
      printf("%f  %f  %f  %f  %f  %f  %f\n",k,a,b,f(a),f(b),p,f(p))
      disp('-------------------------------------------------------------------')
      printf('A raiz está entre %f e %f\n',a,b')
          disp('Pelo método da Bissecção')
      break
      
    else  
      printf("%f  %f  %f  %f  %f  %f  %f\n",k,a,b,f(a),f(b),p,f(p))
      if f(a)*f(p)<0
        b=p;
        
      else
        a=p;
        
      endif
      
      k=k+1;
    endif
    
  endwhile
endif
%%SAÍDA%%
k
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2 %MÉTODO Falsa Posição

%%%INPUT%%%

disp('-------------------------------------------------------------------')
disp('                  MÉTODO EM USO: Falsa Posição                     ')
disp('-------------------------------------------------------------------')
disp('Defina os intervalos do "chute inicial"')

disp('Insira o valor de a:')
disp('Valor de a usado na resolução: 3.111')
disp('')
a = input('');

disp('Insira o valor de b:')
disp('Valor de b usado na resolução: 3.222')
disp('')
b = input('');
disp('-------------------------------------------------------------------')
disp('Me informe a função que deseja aplicar o método da falsa posição:')
disp('Função do exercicio: @(x) sin(x)')
disp('')
f = input('');
disp('-------------------------------------------------------------------')
disp('Defina o erro a ser usado no critério de parada:')
disp('Erro do exercicio: 0.001')
disp('')
erro = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
disp("                        TABELA DE ITERAÇÕES                        ")
disp('-------------------------------------------------------------------')
disp("  k        a         b          f(a)        f(b)       p       f(p)")

if f(a)*f(b) > 0
   error('Ta errado irmão, troca issae sobrinho ~By Bill!!')
else  
  k=0;
  
  while k<100 
        p = (a*abs(f(b))+b*abs(f(a)))/(abs(f(b))+abs(f(a)));
 
    if abs(b-a)<erro %Critério de parada%
      format long
      printf("%f  %f  %f  %f  %f  %f  %f\n",k,a,b,f(a),f(b),p,f(p))
      disp('-------------------------------------------------------------------')
      printf('A raiz está entre %f e %f\n',a,b')
          disp('Pelo método da Falsa Posição')
      break
      
    else  
      printf("%f  %f  %f  %f  %f  %f  %f\n",k,a,b,f(a),f(b),p,f(p))
      if f(a)*f(p)<0
        b=p;
        
      else
        a=p;
        
      endif
      
      k=k+1;
    endif
    
  endwhile
endif
%%SAÍDA%%
k
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 3 %MÉTODO Newton-Raphson

%%%INPUT%%%

disp('-------------------------------------------------------------------')
disp('                  MÉTODO EM USO: Newton-Raphson                    ')
disp('-------------------------------------------------------------------')
disp('Insira o valor de x:')
disp('EX: 1.5')
disp('')
x = input('');
disp('-------------------------------------------------------------------')
disp('Insira a função que será usada:')
disp('EX: @(x) x .^ 3 - 9 * x + 3')
disp('')
f = input('');
disp('-------------------------------------------------------------------')
disp('Me informe a derivada da função que deseja aplicar o método de Newton-Raphson:')
disp('EX: @(x) 3 * x .^ 2 + 9')
disp('')
df = input('');
disp('-------------------------------------------------------------------')
disp('Defina o erro a ser usado no critério de parada:')
disp('EX: 0.003')
disp('')
erro = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
disp('                   Método de Newton-Raphson                        ')
disp('-------------------------------------------------------------------')
disp("  k          x         f(x)          df(x)        p       abs((p-x)/p)")
disp('-------------------------------------------------------------------')
k = 0;
while k<100
  
  p = x - (f(x)/df(x));
  
  if abs(f(p))<erro %Critério de parada%
    break
    % abs(p - x0)<erro
    % abs(p - x0/p)<erro
  endif
  
  printf('%f   %f   %f   %f   %f   %f\n', k,x,f(x),df(x),p,abs((p-x)/p))
  x = p;
  k = k+1;
endwhile
%%SAÍDA%%
disp('-------------------------------------------------------------------')
p
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 4 %MÉTODO da Secante

%%%INPUT%%%

  disp('-------------------------------------------------------------------')
  disp('                      MÉTODO EM USO: Secante                       ')
  disp('-------------------------------------------------------------------')
  disp('Insira o valor para x0:') %0
  disp('EX: 0')
  disp('')
  x0 = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o valor para x1:') %1
  disp('EX: 1')
  disp('')
  x1 = input('');
  disp('-------------------------------------------------------------------')
  disp('Defina a função usando o @(x):')
  disp('EXEMPLO: @(x) x .^ 3 - 9 * x + 3')
  disp('')
  f = input('');
  disp('-------------------------------------------------------------------')
  disp('Defina o erro para o critério de parada:') %0.0005
  disp('EX: 0.0005')
  disp('')
  erro = input('');
  
%%%FIM DO INPUT%%%

  disp('-------------------------------------------------------------------')
  disp('                        Método da Secante                          ')
  disp('-------------------------------------------------------------------') 
  disp("  k          x0         f(x0)        x1        f(x1)     abs(x0-p)")
  
  k=0;
  while k<100
    p=(x0*f(x1)-x1*f(x0))/(f(x1)-f(x0));
    if f(p)==0 || abs(x0-p)<erro %Critério de Parada%
      break
    else 
      
      printf('%f   %f   %f   %f   %f   %f\n', k,x0,f(x0),x1,f(x1),abs(x0-p))
      
      k=k+1;
      x0=x1;
      x1=p;
    endif
  endwhile
  %%SAÍDA%%
  disp('-------------------------------------------------------------------')
  x0
  x1
  k
  disp('-------------------------------------------------------------------')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 5 %MÉTODO Gauss-Seidel

%EXEMPLO
% A = [3 -0.1 -0.2; 0.1 7 -0.3; 0.3 -0.2 10]
% b = [7.85; -19.3; 71.4]
% erro = 0.01
% x = [0; 0; 0]

disp('-------------------------------------------------------------------')
disp('                   MÉTODO EM USO: Gauss-Siedel                     ')
disp('-------------------------------------------------------------------')

%%%INPUT%%%
disp('ATENÇÃO:')
disp('Lembre-se de analisar o CRITÉRIO DE CONVERGÊNCIA ')
disp('Onde os elementos da diagonal principal são maiores em MÓDULO')
disp('Do que a soma dos MÓDULOS dos outros elementos da linha')
disp('EX:')

[5 2 2; 1 3 1; 0 6 8]

disp('Nesse caso, |5|>|2|+|2|,')
disp('Nesse caso, |3|>|1|+|1| e')
disp('Nesse caso, |8|>|6|+|0|')
disp('')
disp('PORTANTO o critério foi satisfeito e isso garante convergência')
disp('Para o método de Gauss-Siedel e Gauss-Jacobi')
disp('')
disp('Caso não seja satifeito, faça PERMUTAÇÃO de linhas da matriz')
disp('')
disp('-------------------------------------------------------------------')

disp('Insira a matriz de coeficientes do sistema linear:')
disp('FORMATO: [valor valor; valor valor]')
disp('EX: [3 -0.1 -0.2; 0.1 7 -0.3; 0.3 -0.2 10]')
disp('')
A = input('');
disp('-------------------------------------------------------------------')
disp('Insira a matriz COLUNA dos termos INDEPENDENTES do sistema linear:')
disp('FORMATO: [valor; valor; valor]')
disp('EX: [7.85; -19.3; 71.4]')
disp('')
b = input('');
disp('-------------------------------------------------------------------')
disp('Insira o ERRO que deseja utilizar:')
disp('O erro está em porcentagem!!')
disp('EX: 0.01')
disp('')
erro = input('');
disp('-------------------------------------------------------------------')
disp('Insira matriz coluna x dos valores iniciais: ')
disp('Para inserir o x "padrão" digite: 0 (zero)')
disp('')
x = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
if x == 0
  x = [0; 0; 0; 0; 0];
endif

n = length(b);
for i = 1:n %Contador de Linha
  I(i,1) = b(i,1)/A(i,i);
  
  for j = 1:n %Contador de Coluna
    if j==i
      D(i,j)=0;
    else
      D(i,j)=-A(i,j)/A(i,i);
    endif
  endfor
endfor

k=0;
while k<100
    x1=x; %Gauss Siedel
    for i=1:n
      x1(i,1)=I(i,1)+D(i,:)*x1;
    endfor
  
  for i=1:n
    ER(i)=abs((x1(i)-x(i))/x1(i))*100;
  endfor
  
  if max(ER)<erro %Critério de parada%
    break
  endif
  k=k+1;
  x=x1;
endwhile
%%SAÍDA%%
x1
k
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 6 %MÉTODO Gauss-Jacobi

%EXEMPLO
% A = [3 -0.1 -0.2; 0.1 7 -0.3; 0.3 -0.2 10]
% b = [7.85; -19.3; 71.4]
% erro = 0.01
% x = [0; 0; 0]

disp('-------------------------------------------------------------------')
disp('                   MÉTODO EM USO: Gauss-Jacobi                     ')
disp('-------------------------------------------------------------------')

%%%INPUT%%%
disp('ATENÇÃO:')
disp('Lembre-se de analisar o CRITÉRIO DE CONVERGÊNCIA ')
disp('Onde os elementos da diagonal principal são maiores em MÓDULO')
disp('Do que a soma dos MÓDULOS dos outros elementos da linha')
disp('EX:')

[5 2 2; 1 3 1; 0 6 8]

disp('Nesse caso, |5|>|2|+|2|,')
disp('Nesse caso, |3|>|1|+|1| e')
disp('Nesse caso, |8|>|6|+|0|')
disp('')
disp('PORTANTO o critério foi satisfeito e isso garante convergência')
disp('Para o método de Gauss-Siedel e Gauss-Jacobi')
disp('')
disp('Caso não seja satifeito, faça PERMUTAÇÃO de linhas da matriz')
disp('')
disp('-------------------------------------------------------------------')

disp('Insira a matriz de coeficientes do sistema linear:')
disp('FORMATO: [valor valor; valor valor]')
disp('EX: [3 -0.1 -0.2; 0.1 7 -0.3; 0.3 -0.2 10]')
disp('')
A = input('');
disp('-------------------------------------------------------------------')
disp('Insira a matriz COLUNA dos termos INDEPENDENTES do sistema linear:')
disp('FORMATO: [valor; valor; valor]')
disp('EX: [7.85; -19.3; 71.4]')
disp('')
b = input('');
disp('-------------------------------------------------------------------')
disp('Insira o ERRO que deseja utilizar:')
disp('O erro está em porcentagem!!')
disp('EX: 0.01')
disp('')
erro = input('');
disp('-------------------------------------------------------------------')
disp('Insira matriz coluna x dos valores iniciais: ')
disp('Para inserir o x "padrão" (matriz nula) digite: 0 (zero)')
disp('')
x = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
if x==0
  x = [0; 0; 0];
endif

n = length(b);
for i = 1:n %Contador de Linha
  I(i,1) = b(i,1)/A(i,i);
  
  for j = 1:n %Contador de Coluna
    if j==i
      D(i,j)=0;
    else
      D(i,j)=-A(i,j)/A(i,i);
    endif
  endfor
endfor

k=0;
while k<100
   x1 = I+D*x; %Gauss Jacobi
  
  for i=1:n
    ER(i)=abs((x1(i)-x(i))/x1(i))*100;
  endfor
  
  if max(ER)<erro %Critério de parada%
    break
  endif
  k=k+1;
  x=x1;
endwhile
%%SAÍDA%%
x1
k
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 7 %Fatoração LU

disp('-------------------------------------------------------------------')
disp('                   MÉTODO EM USO: Fatoração LU                     ')
disp('-------------------------------------------------------------------')


%%%INPUT%%%

disp('Coloque a matriz A para a fatoração LU:')
disp('EX: [3 -4 1; 1 2 2; 4 0 -3]')
disp('')
A = input('');

disp('Coloque a matriz b dos termos independentes:')
disp('EX: [9; 3; -2]')
disp('')
b = input('');

%%%FIM DO INPUT%%%

disp('-------------------------------------------------------------------')
[L,U,P] = lu(A);
y = linsolve(L,P*b)
x = inv(U)*y
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 8 %MÉTODOS de EDO
opcao2 = menu('MÉTODOS de EDO', '1)Método de Euler', '2)Método do Ponto Médio', '3)Métodos de Runge-Kutta', '4)EDOs de ordem superior')
    
  switch(opcao2)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1 %Método de Euler
  
  disp('-------------------------------------------------------------------')
  disp('                       MÉTODO EM USO: Euler                        ')
  disp('-------------------------------------------------------------------')

  %%%INPUT%%%
  
  disp('Insira a EDO:')
  disp('EX: @(x,y) -2*(x^3)+12*(x^2)-20*x+8.5+0*y')
  disp('')
  edo = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o valor inicial de X, chamado de x0:')
  disp('EX: 0')
  disp('')
  x0 = input('');
  
  disp('Insira o valor inicial de Y, chamado de y0:')
  disp('EX: 1')
  disp('')
  y0 = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o passo que será utilizado, chamado por h:')
  disp('EX: 0.5')
  disp('')
  h = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira até qual x você deseja chegar, chamado de xn:')
  disp('EX: 4')
  disp('')
  xn = input('');
  disp('-------------------------------------------------------------------')
  
  %%%FIM DO INPUT%%%
  disp("                         RESOLUÇÃO DA EDO                          ")
  disp('-------------------------------------------------------------------')
  disp("    i        x(i)        y(i)       k1")
  
    y(1) = y0;
  x = x0:h:xn;
  n = length(x);
  i = 1;

  while i < (n+1)
    K1 = edo(x(i),y(i));
    y(i+1) = y(i) + K1*h;
    
    printf('%f   %f   %f   %f\n',i, x(i),y(i),K1)
    i = i+1;
  endwhile
  %%SAÍDA%%
  
  disp('-------------------------------------------------------------------')
  disp('A resolução da EDO dentro do "intervalo" estupulado são os valores ')
  disp('de X e Y na tabela')
  disp('-------------------------------------------------------------------')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 2 %Método do Ponto Médio
  
  disp('-------------------------------------------------------------------')
  disp('                   MÉTODO EM USO: Ponto Médio                      ')
  disp('-------------------------------------------------------------------')

  %%%%%INPUT%%%%%
  
  disp('-------------------------------------------------------------------')
  disp('Insira a EDO:')
  disp('')
  edo = input('');
  disp('-------------------------------------------------------------------')
  disp('Insira o valor inicial de X, chamado de x0:')
  disp('')
  x0 = input('');
  
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
  xn = input('');
  disp('-------------------------------------------------------------------')
  
  %%%FIM DO INPUT%%%
  disp("  x(i)     y(i)")
  
  n = (xn-x0)/h;

  x(1) = x0;
  y(1) = y0;

  for i = 1:n
    
    x(i+1) = x(i)+h;
    K1 = edo(x(i),y(i));
    
    ym = y(i) + K1 * (h/2); %retira para euler%
    K2 = edo(x(1)+(h/2),ym); %retira para euler%
    
    %y(i+1) = y(i)+K1*h %troca com a de baixo para euler%
    y(i+1) = y(i)+K2*h %retira para euler%
    printf('%f   %f\n', x(i),y(i))
  endfor
  %%SAÍDA%%
  
  disp('-------------------------------------------------------------------')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 3 %Métodos de Runge-Kutta
    opcao3 = menu('MÉTODOS DE RUNGE-KUTTA', '1)Runge-Kutta 1ºORDEM', '2)Runge-Kutta 2ºORDEM', '3)Método de Ralson', '4)Runge-Kutta 3ºORDEM','5)Runge-Kutta 4ºORDEM')
    
    switch(opcao3)
    
    case 1 %Runge-Kutta 1ºORDEM%
      disp('O método de Runge-Kutta de 1ºOrdem nada mais é que o Método de Euler')
      disp('')
      disp('Deseja rodá-lo?')
      disp('1)SIM')
      disp('2)NÃO')
      rk1 = input('');
      
      if rk1 == 1
        disp('-------------------------------------------------------------------')
        disp('            MÉTODO EM USO: Euler ou Runge-Kutta 1ºORDEM            ')
        disp('-------------------------------------------------------------------')

        %%%INPUT%%%
        
        disp('Insira a EDO:')
        disp('EX: @(x,y) -2*(x^3)+12*(x^2)-20*x+8.5+0*y')
        disp('')
        edo = input('');
        disp('-------------------------------------------------------------------')
        disp('Insira o valor inicial de X, chamado de x0:')
        disp('EX: 0')
        disp('')
        x0 = input('');
        
        disp('Insira o valor inicial de Y, chamado de y0:')
        disp('EX: 1')
        disp('')
        y0 = input('');
        disp('-------------------------------------------------------------------')
        disp('Insira o passo que será utilizado, chamado por h:')
        disp('EX: 0.5')
        disp('')
        h = input('');
        disp('-------------------------------------------------------------------')
        disp('Insira até qual x você deseja chegar, chamado de xn:')
        disp('EX: 4')
        disp('')
        xn = input('');
        disp('-------------------------------------------------------------------')
        
        %%%FIM DO INPUT%%%
        disp("                         RESOLUÇÃO DA EDO                          ")
        disp('-------------------------------------------------------------------')
        disp("    i        x(i)        y(i)       k1")
        
          y(1) = y0;
        x = x0:h:xn;
        n = length(x);
        i = 1;

        while i < (n+1)
          K1 = edo(x(i),y(i));
          y(i+1) = y(i) + K1*h;
          
          printf('%f   %f   %f   %f\n',i, x(i),y(i),K1)
          i = i+1;
        endwhile
        %%SAÍDA%%
        
        disp('-------------------------------------------------------------------')
        disp('A resolução da EDO dentro do "intervalo" estupulado são os valores ')
        disp('de X e Y na tabela')
        disp('-------------------------------------------------------------------')
      endif
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 2 %Runge-Kutta 2ºORDEM%
      disp('O método de Runge-Kutta de 2ºOrdem nada mais é que o Método do Ponto Médio')
      disp('')
      disp('Deseja rodá-lo?')
      disp('1)SIM')
      disp('2)NÃO')
      rk2 = input('');
      
      if rk2 == 1
        disp('-------------------------------------------------------------------')
        disp('         MÉTODO EM USO: Ponto Médio ou Runge-Kutta 2ºORDEM         ')
        disp('-------------------------------------------------------------------')

        %%%%%INPUT%%%%%
        
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
        disp("  x(i)     y(i)")
        
        n = (Xn-X0)/h;

        X(1) = X0;
        y(1) = y0;

        for i = 1:n
          
          X(i+1) = X(i)+h;
          K1 = edo(X(i),y(i));
          
          ym = y(i) + K1 * (h/2); %retira para euler%
          K2 = edo(X(1)+(h/2),ym); %retira para euler%
          
          %y(i+1) = y(i)+K1*h %troca com a de baixo para euler%
          y(i+1) = y(i)+K2*h;%% %retira para euler%
          printf('%f   %f\n', X(i),y(i))
        endfor
        %%SAÍDA%%
        
        disp('-------------------------------------------------------------------')
      endif
    case 3 %Método de Ralson%
      
      disp('-------------------------------------------------------------------')
      disp('                   EDO pelo método de Ralson                       ')
      disp('-------------------------------------------------------------------')

      disp('Insira a EDO:')
      disp('EX: @(x,y) 0.2*10*((100*(20-25)*(20-x))/100-2.5*x)+0*y')
      disp('')
      edo1 = input('');

      disp('-------------------------------------------------------------------')
      disp('Insira o valor inicial de X, chamado de x0:')
      disp('EX: 0')
      disp('')
      x0 = input('');

      disp('Insira o valor inicial de Y, chamado de y0:')
      disp('EX: 0.7')
      disp('')
      y0 = input('');

      %disp('Insira o valor inicial de Z, chamado de z0:')
      %disp('Exercício: NÃO TINHA')
      %disp('')
      %z0 = input('');

      disp('-------------------------------------------------------------------')
      disp('Insira o passo que será utilizado, chamado por h:')
      disp('EX: 0.1')
      disp('')
      h = input('');

      disp('-------------------------------------------------------------------')
      disp('Insira até qual x você deseja chegar, chamado de xn:')
      disp('EX: 1')
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
        disp('Vale destacar que o método de Ralson é um dos métodos de Runge-Kutta')
        disp('de segunda ordem, sendo o método com o MENOR erro de truncamento.')
        disp('-------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 4 %Runge-Kutta 3ºORDEM%
      
      disp('-------------------------------------------------------------------')
      disp('         MÉTODO EM USO: Runge-Kutta 3ºORDEM         ')
      disp('-------------------------------------------------------------------')

      %%%%%INPUT%%%%%
  
      disp('Insira a EDO:')
      disp('')
      edo = input('');
      disp('-------------------------------------------------------------------')
      disp('Insira o valor inicial de X, chamado de x0:')
      disp('')
      x0 = input('');
        
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
      xn = input('');
      disp('-------------------------------------------------------------------')
        
      %%%FIM DO INPUT%%%
      disp("  x(i)     y(i)")
      
      n=(xn - x0)/h;
      x(1) = x0;
      y(1) = y0;
      
      for i = 1:n
        x(i+1) = x(i) + h;
        K1 = edo(x(i),y(i));
        K2 = edo(x(i) + (h/2), y(i) + K1 * (h/2));
        K3 = edo(x(i) + h, y(i) - (K1 * h) + (2 * K2 * h));
        y(i+1) = y(i) + (K1 + 4*K2 + K3) * (h/6); %%
        printf('%f   %f\n', x(i),y(i))
      endfor
      %%SAÍDA%%
      
      disp('-------------------------------------------------------------------')
    case 5 %Runge-Kutta 4ºORDEM%
      
      disp('-------------------------------------------------------------------')
      disp('         MÉTODO EM USO: Runge-Kutta 4ºORDEM         ')
      disp('-------------------------------------------------------------------')

      %%%%%INPUT%%%%%
  
      disp('Insira a EDO:')
      disp('')
      edo = input('');
      disp('-------------------------------------------------------------------')
      disp('Insira o valor inicial de X, chamado de x0:')
      disp('')
      x0 = input('');
        
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
      xn = input('');
      disp('-------------------------------------------------------------------')
        
      %%%FIM DO INPUT%%%
      disp("  x(i)     y(i)")
      
      n=(xn - x0)/h;
      x(1) = x0;
      y(1) = y0;
      
      for i = 1:n
        x(i+1) = x(i) + h;
        K1 = edo(x(i),y(i));
        K2 = edo(x(i) + (h/2), y(i) + K1 * (h/2));
        K3 = edo(x(i) + (h/2), y(i) + K2 * (h/2));
        K4 = edo(x(i) + h, y(i) + K3 *h);
        y(i+1) = y(i) + (K1 + 2*K2 + 2*K3 + K4) * (h/6); %%
        printf('%f   %f\n', x(i),y(i))
      endfor
      %%SAÍDA%%
      
      disp('-------------------------------------------------------------------')
    endswitch %switch do Runge-Kutta%
    
    %disp('Deseja usar o programa novamente?')
    %disp('1)SIM')
    %disp('2)NÃO')
    %disp('')
    %repe = input('');

    %if repe == 1
    %  oneforall
    %else
    %  disp('-------------------------------------------------------------------')
    %  disp('Até Mais, yeeey!!')
    %  disp('Qualquer coisa estou aqui...')
    %  disp('U.u   zZZZZzzzzzZzzzzZZz')
    %endif

  case 4 %Ordem Superior%
  
    disp('-------------------------------------------------------------------')
    disp('               MÉTODO EM USO: EDOs de ordem superior               ')
    disp('-------------------------------------------------------------------')

    %%%%%INPUT%%%%%
  
    disp('Insira a EDO Y:')
    disp('')
    edoy = input('');
    disp('Insira a EDO Z:')
    disp('')
    edoz = input('');
    disp('-------------------------------------------------------------------')
    disp('Insira o valor inicial de X, chamado de x0:')
    disp('')
    x0 = input('');
    
    disp('Insira o valor inicial de Y, chamado de y0:')
    disp('')
    y0 = input('');
    
    disp('Insira o valor inicial de Z, chamado de y0:')
    disp('')
    z0 = input('');
    disp('-------------------------------------------------------------------')
    disp('Insira o passo que será utilizado, chamado por h:')
    disp('')
    h = input('');
    disp('-------------------------------------------------------------------')
    disp('Insira até qual x você deseja chegar, chamado de xn:')
    disp('')
    xn = input('');
    disp('-------------------------------------------------------------------')
        
    %%%FIM DO INPUT%%%
    disp("  x(i)     y(i)     z(i)")
    
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
      printf('%f   %f   %f\n', x(i), y(i), z(i))
    endfor
    %%SAÍDAS%%
    
    disp('-------------------------------------------------------------------')
  
disp('-------------------------------------------------------------------')
endswitch %Switch de EDOs%

%disp('Deseja usar o programa novamente?')
%disp('1)SIM')
%disp('2)NÃO')
%disp('')
%repe = input('');

%if repe == 1
% oneforall
%else
%  disp('-------------------------------------------------------------------')
%  disp('Até Mais, yeeey!!')
%  disp('Qualquer coisa estou aqui...')
%  disp('U.u   zZZZZzzzzzZzzzzZZz')
%endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 9 %Interpolação Polinomial

  opcao3 = menu('Interpolação Polinomial:', '1)MÉTODO: Lagrange', '2)MÉTODO: Newton/DIF.DIV.')
  
  switch(opcao3)
  
    case 1 %Lagrange%
    
      disp('--------------------------------------------')
      disp('           Método de Lagrange               ')
      disp('--------------------------------------------')

      %INPUT%

      disp('Insira a matriz dos valores de X')
      disp('EX: [-1 0 1]')
      x = input('');
      disp('--------------------------------------------')
      disp('Insira a matriz dos valores de Y')
      disp('EX: [6 1 0]')
      y = input('');
      disp('--------------------------------------------')
      disp('Insira o valor de X que será interpolado')
      disp('EX: 0.9')
      xint = input('');

      %FIM DO INPUT%

      n = length(x);

      for i = 1:n
         L(i)= 1;
         
         for j = 1:n
           if i != j
             L(i) = L(i)*((xint-x(j))/(x(i)-x(j)));
             %L do Exemplo: -0.045    0.19    0.855%
           endif
         endfor
      endfor

      yint = 0; %Resposta do exemplo: -0.08%

      for i = 1:n
        yint = yint + y(i)*L(i);
      endfor
      disp('--------------------------------------------')
      printf('O valor de Y para o X interpolado é igual a: \n %f \n',yint)
      disp('--------------------------------------------')
      disp('          VALOR DAS Ls CALCULADAS           ')
      disp('--------------------------------------------')
      L
      disp('--------------------------------------------')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 %Newton/DIF.DIV.%
    
      disp('--------------------------------------------')
      disp('       Método de Newton INTERPOLAÇÃO        ')
      disp('--------------------------------------------')

      %INPUT%

      disp('Insira a matriz dos valores de X')
      disp('EX: [-1 0 1]')
      x = input('');
      disp('--------------------------------------------')
      disp('Insira a matriz dos valores de Y')
      disp('EX: [6 1 0]')
      y = input('');
      disp('--------------------------------------------')
      disp('Insira o valor de X que será interpolado')
      disp('EX: 0.9')
      xint = input('');

      %FIM DO INPUT%

      n = length(x);

      for i = 1:n
        DD(i,1) = y(i);
      endfor

      for j = 2:n
        for i = 1:(n-j+1)
          DD(i,j) = (DD(i+1,j-1) - DD(i,j-1))/(x(i+j-1) - x(i));
          %DD do exemplo [6 -5 2; 1 -1 0; 0 0 0]%
        endfor
      endfor

      for j = 1:n
        a(j) = DD(1,j);
      endfor

      xtemp = 1;
      yint = a(1);

      for k = 1:n-1
        xtemp = xtemp * (xint - x(k));
        yint = yint + a(k+1) * xtemp;  
      endfor

      disp('--------------------------------------------')
      printf('O valor de Y para o X interpolado é igual a: \n %f \n',yint)
      disp('--------------------------------------------')
      disp('          VALOR DAS Ls CALCULADAS           ')
      disp('--------------------------------------------')
      DD
      disp('--------------------------------------------')      
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  endswitch

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 10 %Integração Numérica
  
disp('-----------------------------------------------------------------')  
disp('Opção corrompida...')
disp('ou não implementada.....por preguiça....')
disp('-----------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 11 %Ajuste de Curvas

%   AJUSTE DE CURVAS   %
%  Versão 2.9kb  %

linear = menu('Ajuste de curvas:', '1)Ajuste Linear','2)Ajuste Polinomial', '3)Ajuste não linear','4)Ajuste Customizado', '5)Gráfico de Dispersão')
switch(linear)
  case 1 %  1)Ajuste Linear
  
    disp('-----------------------------------------------------------------')
    disp('                          AJUSTE LINEAR                          ')
    disp('-----------------------------------------------------------------')
    disp('Insira os valores de X, dos pares ordenados em formato de MATRIZ')
    disp('EX: [183 173 168 188 158 163 193 163 178]') %[0.5 0.75 1 1.5 2 2.5 3]
    disp('')
    x = input('');  % Matriz de X

    disp('')

    disp('Insira os valores de Y, dos pares ordenados em formato de MATRIZ')
    disp('EX: [79 69 70 81 61 63 79 71 73]') %[-2.8 -0.6 1 3.2 4.8 6 7]
    disp('')
    y = input('');  % Matriz de Y

    disp('')

    disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
    disp('EX: 175')
    disp('')
    xe = input('');  % Valor de X teste "estimado"

    disp('')

    %  CASO PRECISE "INVERTER" X COM Y, BASTA COLOCAR PRIMEIRO A MATRIZ DE Y E DEPOIS DE X, COLOCANDO UM VALOR ESTIMADO DE Y  %

    %disp('Insira um valor de Y para ser utilizado como teste para o ajuste:')
    %disp('EX: 80')
    %disp('')
    %ye = input('');  % Valor de Y teste "estimado"

    %  FIM DO INPUT  %

    m = length(x);

    %  SOMATÓRIOS  %
    sx = sum(x);
    sy = sum(y);
    sxQUAD = sum(x.^2);
    syQUAD = sum(y.^2);
    sxy = sum(x.*y);
    sx2 = (sum(x))^2;
    sy2 = (sum(y))^2;
    %              %

    disp('-----------------------------------------------------------------')
    disp('Deseja visualisar o valor dos somátórios?')
    disp('1)Sim')
    disp('2)Nao')
    disp('')
    versomats = input('');

    if versomats == 1
      disp('-----------------------------------------------------------------')
      printf('Somatório de x: %f\n',sx)
      printf('Somatório de y: %f\n',sy)
      printf('Somatório de x^2: %f\n',sxQUAD)
      printf('Somatório de x*y: %f\n',sxy)
    endif


    %  Cálculo dos Alfas, em somatório  %

    ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);

    ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);

    %  R^2  %
    disp('-----------------------------------------------------------------')
    disp('A equação da reta encontrada com os valores de Alfa 0 e Alfa 1 é:')
    disp('-----------------------------------------------------------------')
    printf('y = %f + (%f) * x\n',ALFA0, ALFA1)

    disp('')

    printf('Onde "Alfa 0" é igual a: %f\n',ALFA0) % Alfa 0 do exemplo: -20.078
    disp('e')
    printf('Onde "Alfa 1" é igual a: %f\n',ALFA1) % Alfa 1 do exemplo: 0.5276


    f0 = @(x) ALFA0 + ALFA1 * x;    % EQUAÇÃO DA RETA COM PARES DADOS
    f1 = @(xe) ALFA0 + ALFA1 * xe;  % EQUAÇÃO DA RETA COM X DE TESTE

    %  COEFICIENTE DE DETERMINAÇÃO  %
    ymedio = sy/m;

    St = sum((y - ymedio).^2); % Média

    Sr = sum((y - f0(x)).^2);  % Curva ajustada

    RQUAD = (St-Sr)/St; % R^2 do exemplo: 0.8366

    RQuadComputacional = (((m*sxy)-(sx*sy))/(sqrt((m*sxQUAD)-sx2)*(sqrt((m*syQUAD)-sy2))))^2; % R^2 computacional funciona para o linear 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('-----------------------------------------------------------------')
    disp('                 Coeficiente de determinação:                    ')
    disp('-----------------------------------------------------------------')

    printf('Valor do Y médio: %f\n',ymedio)
    printf('Valor do St: %f\n',St)
    printf('Valor do Sr: %f\n',Sr)

    disp('-----------------------------------------------------------------')
    disp('                      Valor do R qudradrado:                     ')
    disp('-----------------------------------------------------------------')

    printf('Valor do R^2: %f\n',RQUAD)

    printf('Valor do R^2 Computacional: %f\n',RQuadComputacional)

    disp('-----------------------------------------------------------------')
    disp('                  Valor de Y para o X estimado:                ')
    disp('-----------------------------------------------------------------')

    printf('O valor de X estimado foi: %f\n',xe)
    printf('O valor de Y para o X estimado foi: %f\n',f1(xe)')

    disp('-----------------------------------------------------------------')
    disp('Deseja visualizar o diagrama de dispersão?')
    disp('1)Sim')
    disp('2)Não')
    disp('')
    diagrama = input('')

    if diagrama == 1
      scatter(x,y)
    endif

    disp('-----------------------------------------------------------------')
    disp('                     FIM DO AJUSTE LINEAR                        ')
    disp('-----------------------------------------------------------------')
  
  case 2 %  2)Ajuste Polinomial
    disp('-----------------------------------------------------------------')
    disp('                      AJUSTE POLINOMIAL                          ')
    disp('-----------------------------------------------------------------')
    disp('Insira o valor em X dos nós (pares ordenados) em formato de MATRIZ')
    disp('EX: [300 350 400]')
    disp('')
    X = input('');

    disp('')
    
    disp('Insira o valor em Y dos nós (pares ordenados) em formato de MATRIZ')
    disp('EX: [1.132 0.978 0.854]')
    disp('')
    Y = input('');
    
    disp('')
    
    disp('-----------------------------------------------------------------')
    disp('A ordem do polinômio será calculada automaticamente, com base')
    disp('na quantidade de pares ordenados de nós que foram inseridos.')
    n = (length(X)-1);
    
    disp('')
    
    printf('O polinômio atual é de %fº ordem!\n',n)
    disp('-----------------------------------------------------------------')
    
    disp('Deseja alterar a ordem do polinômio?')
    disp('1)Sim')
    disp('2)Não')
    alterpoli = input('');
    
    if alterpoli == 1
      disp('Insira a nova ordem do polinômio que será ajustado:')
      disp('')
      n = input('');
      if n < (length(X)-1)
        error('A ordem do polinômio não pode ser menor que: %f',(length(X)-1))
      endif
    endif
    
    p = polyfit(X,Y,n);
    
    disp('-----------------------------------------------------------------')
    disp('O polinômio obtido foi:')
    polyout(p)
    disp('-----------------------------------------------------------------')
    disp('                    FIM DO AJUSTE POLINOMIAL                     ')
    disp('-----------------------------------------------------------------')
    
  case 3 %  3)Ajuste não linear
    disp('-----------------------------------------------------------------')
    disp('                       AJUSTE NÃO LINEAR                         ')
    disp('-----------------------------------------------------------------')
    disp('Insira os valores de X e Y para visualizar o diagrama de dispersão:')
    disp('INSIRA X, no formato de MATRIZ:')
    disp('EX Exponencial: [-1 -0.7 -0.4 -0.1 0.2 0.5 0.8 1]')
    disp('')
    
    x = input('');
    disp('INSIRA Y, no formato de MATRIZ:')
    disp('EX Exponencial: [36.547 17.264 8.155 3.852 1.820 0.860 0.406 0.246]')
    disp('')
    
    y = input('');
    
    m = length(x);
    
    scatter(x,y)
    
    disp('-----------------------------------------------------------------')
    disp('                ANALISE O GRÁFICO DE DISPERSÃO                   ')
    disp('-----------------------------------------------------------------')
    disp('Com base no gráfico de disperção, escolha qual função descreve   ')
    disp('melhor a dispersão dos valores inseridos.')
    disp('-----------------------------------------------------------------')
    
    nlinear = menu('Não lineares:', '1)Função Exponencial', '2)Função Logarítmica', '3)Função Potencial', '4)Função Hiperbólica')
    switch(nlinear)
      case 1 %  1)Função Exponencial
     
        %  Exponencial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE NÃO LINEAR: EXPONENCIAL                   ')
        disp('-----------------------------------------------------------------')
        
        Y = log(y);
        X = x;
        
        sx = sum(X);
        sy = sum(Y);
        sxQUAD = sum(X.^2);
        syQUAD = sum(Y.^2);
        sxy = sum(X.*Y);
        sx2 = (sum(X))^2;
        sy2 = (sum(Y))^2;
        
        disp('Deseja visualisar o valor dos somátórios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somatório de X: %f\n',sx)
          printf('Somatório de Y: %f\n',sy)
          printf('Somatório de X^2: %f\n',sxQUAD)
          printf('Somatório de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = exp(ALFA0);
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A função exponencial ajustada para esses valores é:')
        printf('Y = %f * exp(%f * X)\n',a,b)
       
        disp('')

        printf('Onde "a" é igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" é igual a: %f\n',b) % b do exemplo:
       
        fexp0 = @(X) a * exp(b * X); % EQUAÇÃO EXPONENCIAL COM PARES DADOS
        fexp1 = @(Xe) a * exp(b * Xe); % EQUAÇÃO EXPONENCIAL COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINAÇÃO  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % Média

        Sr = sum((Y - fexp0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gráfico de Dispesão                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        %disp('-----------------------------------------------------------------')
        %disp('                 Coeficiente de determinação:                    ')
        %disp('-----------------------------------------------------------------')

        %printf('Valor do Y médio: %f\n',ymedio)
        %printf('Valor do St: %f\n',St)
        %printf('Valor do Sr: %f\n',Sr)

        %disp('-----------------------------------------------------------------')
        %disp('                      Valor do R qudradrado:                     ')
        %disp('-----------------------------------------------------------------')

        %printf('Valor do R^2: %f\n',RQUAD)
        
        disp('Deseja testar a função encontrada?')
        disp('1)Sim')
        disp('2)Não')
        disp('')
        testarfun = input('')
        
        if testarfun == 1
          disp('-----------------------------------------------------------------')
          disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
          disp('')
          Xe = input('');  % Valor de X teste "estimado"
          disp('-----------------------------------------------------------------')
          disp('                  Valor de Y para o X estimado:                ')
          disp('-----------------------------------------------------------------')

          printf('O valor de X estimado foi: %f\n',Xe)
          printf('O valor de Y para o X estimado foi: %f\n',fexp0(Xe)')
        endif
        
        disp('-----------------------------------------------------------------')
        disp('             FIM DO AJUSTE NÃO LINEAR: EXPONENCIAL               ')
        disp('-----------------------------------------------------------------')
        
      case 2 %  2)Função Logarítmica %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE NÃO LINEAR: LOGARÍTMICA                   ')
        disp('-----------------------------------------------------------------')
        
        Y = y;
        X = log(x);
        
        sx = sum(X);
        sy = sum(Y);
        sxQUAD = sum(X.^2);
        syQUAD = sum(Y.^2);
        sxy = sum(X.*Y);
        sx2 = (sum(X))^2;
        sy2 = (sum(Y))^2;
        
        disp('Deseja visualisar o valor dos somátórios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somatório de X: %f\n',sx)
          printf('Somatório de Y: %f\n',sy)
          printf('Somatório de X^2: %f\n',sxQUAD)
          printf('Somatório de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = ALFA1;
        b = exp(ALFA0/a);
        
        disp('-----------------------------------------------------------------')
        disp('A função logarítmica ajustada para esses valores é:')
        printf('Y = %f * log(%f * X)\n',a,b)
       
        disp('')

        printf('Onde "a" é igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" é igual a: %f\n',b) % b do exemplo:
       
        flogarit0 = @(X) a * log(b * X); % EQUAÇÃO LOGARÍTMICA COM PARES DADOS
        flogarit1 = @(Xe) a * log(b * Xe); % EQUAÇÃO LOGARÍTMICA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINAÇÃO  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % Média

        Sr = sum((Y - flogarit0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gráfico de Dispesão                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a função encontrada?')
        disp('1)Sim')
        disp('2)Não')
        disp('')
        testarfun = input('')
        
        if testarfun == 1
          disp('-----------------------------------------------------------------')
          disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
          disp('')
          Xe = input('');  % Valor de X teste "estimado"
          disp('-----------------------------------------------------------------')
          disp('                  Valor de Y para o X estimado:                ')
          disp('-----------------------------------------------------------------')

          printf('O valor de X estimado foi: %f\n',Xe)
          printf('O valor de Y para o X estimado foi: %f\n',flogarit0(Xe)')
        endif
        
        disp('-----------------------------------------------------------------')
        disp('             FIM DO AJUSTE NÃO LINEAR: LOGARÍTMICA               ')
        disp('-----------------------------------------------------------------')
      
      case 3 %  3)Função Potencial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE NÃO LINEAR: POTENCIAL                   ')
        disp('-----------------------------------------------------------------')
        
        Y = log(y);
        X = log(x);
        
        sx = sum(X);
        sy = sum(Y);
        sxQUAD = sum(X.^2);
        syQUAD = sum(Y.^2);
        sxy = sum(X.*Y);
        sx2 = (sum(X))^2;
        sy2 = (sum(Y))^2;
        
        disp('Deseja visualisar o valor dos somátórios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somatório de X: %f\n',sx)
          printf('Somatório de Y: %f\n',sy)
          printf('Somatório de X^2: %f\n',sxQUAD)
          printf('Somatório de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = exp(ALFA0);
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A função potencial ajustada para esses valores é:')
        printf('Y = %f * (X^%f)\n',a,b)
       
        disp('')

        printf('Onde "a" é igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" é igual a: %f\n',b) % b do exemplo:
       
        fpoten0 = @(X) a * (X^b); % EQUAÇÃO POTENCIAL COM PARES DADOS
        fpoten1 = @(Xe) a * (Xe^b); % EQUAÇÃO POTENCIAL COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINAÇÃO  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % Média

        Sr = sum((Y - fpoten0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gráfico de Dispesão                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a função encontrada?')
        disp('1)Sim')
        disp('2)Não')
        disp('')
        testarfun = input('')
        
        if testarfun == 1
          disp('-----------------------------------------------------------------')
          disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
          disp('')
          Xe = input('');  % Valor de X teste "estimado"
          disp('-----------------------------------------------------------------')
          disp('                  Valor de Y para o X estimado:                ')
          disp('-----------------------------------------------------------------')

          printf('O valor de X estimado foi: %f\n',Xe)
          printf('O valor de Y para o X estimado foi: %f\n',fpoten0(Xe)')
        endif
        
        disp('-----------------------------------------------------------------')
        disp('              FIM DO AJUSTE NÃO LINEAR: POTENCIAL                ')
        disp('-----------------------------------------------------------------')
      
      case 4 %  4)Função Hiperbólica  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE NÃO LINEAR: HIPERBÓLICA                   ')
        disp('-----------------------------------------------------------------')
        
        Y = y;
        X = 1/x;
        
        sx = sum(X);
        sy = sum(Y);
        sxQUAD = sum(X.^2);
        syQUAD = sum(Y.^2);
        sxy = sum(X.*Y);
        sx2 = (sum(X))^2;
        sy2 = (sum(Y))^2;
        
        disp('Deseja visualisar o valor dos somátórios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somatório de X: %f\n',sx)
          printf('Somatório de Y: %f\n',sy)
          printf('Somatório de X^2: %f\n',sxQUAD)
          printf('Somatório de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = ALFA0;
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A função hiperbólica ajustada para esses valores é:')
        printf('Y = %f + (%f/x)\n',a,b)
       
        disp('')

        printf('Onde "a" é igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" é igual a: %f\n',b) % b do exemplo:
       
        fhiperbo0 = @(X) a + (b/X); % EQUAÇÃO HIPERBÓLICA COM PARES DADOS
        fhiperbo1 = @(Xe) a + (b/Xe); % EQUAÇÃO HIPERBÓLICA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINAÇÃO  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % Média

        Sr = sum((Y - fhiperbo0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gráfico de Dispesão                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a função encontrada?')
        disp('1)Sim')
        disp('2)Não')
        disp('')
        testarfun = input('')
        
        if testarfun == 1
          disp('-----------------------------------------------------------------')
          disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
          disp('')
          Xe = input('');  % Valor de X teste "estimado"
          disp('-----------------------------------------------------------------')
          disp('                  Valor de Y para o X estimado:                ')
          disp('-----------------------------------------------------------------')

          printf('O valor de X estimado foi: %f\n',Xe)
          printf('O valor de Y para o X estimado foi: %f\n',fhiberpo0(Xe)')
        endif
        
        disp('-----------------------------------------------------------------')
        disp('             FIM DO AJUSTE NÃO LINEAR: HIPERBÓLICA               ')
        disp('-----------------------------------------------------------------')
    endswitch
  
  case 4 %  4)Ajuste Customizado  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE NÃO LINEAR: CUSTOMIZADO                   ')
        disp('-----------------------------------------------------------------')
        disp('')
        disp('OBSERVAÇÕES:')
        disp('OBS 1: Faça a análise da sua função não linear e linearize-a.    ')
        disp('OBS 2: Identifique quem é o termo INDEPENDENTE e o DEPENDENTE.   ')
        disp('OBS 3: O termo DEPENDENTE estará em função do INDEPENDENTE.      ')
        disp('Por exemplo: Y(x) "Y em função de X"                             ')
        disp('OBS 4: Idenfique quem será Y e X.                                ')
        disp('OBS 5: Identifique quem será (se existir) a e b.                 ')
        disp('OBS 6: EDITE esse trecho do código, para sua linearização customizada')
        disp('Basta pesquisar no editor: CUSTOMIZEAQUI')
        disp('')
        disp('-----------------------------------------------------------------')
        disp('                       EDITE O PROGRAMA                          ')
        disp('-----------------------------------------------------------------')
        
        % Tire do modo COMENTÁRIO, as linhas com %---%
        
        %disp('Defina Y comparando com o y da equação que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('Y = log(y)')
        %disp('')
        
        %---%Y = %EDITE AQUI;
        
        %disp('Defina X comparando com o x da equação que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('X = log(x)')
        %disp('')
        
        %---%X = %EDITE AQUI;
        
        %%%%%%%%%%%%%%%%%%%%
        
        %---%sx = sum(X);
        %---%sy = sum(Y);
        %---%sxQUAD = sum(X.^2);
        %---%syQUAD = sum(Y.^2);
        %---%sxy = sum(X.*Y);
        %---%sx2 = (sum(X))^2;
        %---%sy2 = (sum(Y))^2;
        
        %---%disp('Deseja visualisar o valor dos somátórios?')
        %---%disp('1)Sim')
        %---%disp('2)Nao')
        %---%disp('')
        %---%versomats = input('');

        %---%if versomats == 1
          %---%disp('-----------------------------------------------------------------')
          %---%printf('Somatório de X: %f\n',sx)
          %---%printf('Somatório de Y: %f\n',sy)
          %---%printf('Somatório de X^2: %f\n',sxQUAD)
          %---%printf('Somatório de X*Y: %f\n',sxy)
        %---%endif
        
        %---%ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        %---%ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        %---%disp('-----------------------------------------------------------------')
        %disp('Defina a comparando com o alfa0 da equação que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('a = exp(ALFA0)')
        %disp('')
        
        %---%a = %EDITE AQUI;
        
        %disp('Defina b comparando com o alfa1 da equação que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('b = ALFA1')
        %disp('')
        
        %---%b = %EDITE AQUI;
        
        %---%disp('-----------------------------------------------------------------')
        %---%disp('A função customuzada ajustada para esses valores é:')
        %---%printf('Y = EDITE AQUI %f e %f\n',a,b)
       
        %---%disp('')

        %---%printf('Onde "a" é igual a: %f\n',a) % a do exemplo:
        %---%disp('e')
        %---%printf('Onde "b" é igual a: %f\n',b) % b do exemplo:
       
        %fcustom0 = @(X) CUSTOMIZADA; % EQUAÇÃO CUSTOMIZADA COM PARES DADOS
        %fcustom1 = @(Xe) CUSTOMIZADA; % EQUAÇÃO CUSTOMIZADA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINAÇÃO  %
        %---%ymedio = sy/m;

        %---%St = sum((Y - ymedio).^2); % Média

        %---%Sr = sum((Y - fcustom0(X)).^2);  % Curva ajustada

        %---%RQUAD = (St-Sr)/St; 
         
        %---%disp('-----------------------------------------------------------------')
        %---%disp('                   Novo Gráfico de Dispesão                      ')
        %---%disp('-----------------------------------------------------------------')
        %---%scatter(X,Y)
        
        %---%disp('Deseja testar a função encontrada?')
        %---%disp('1)Sim')
        %---%disp('2)Não')
        %---%disp('')
        %---%testarfun = input('')
        
        %---%if testarfun == 1
          %---%disp('-----------------------------------------------------------------')
          %---%disp('Insira um valor de X para ser utilizado como teste para o ajuste:')
          %---%disp('')
          %---%Xe = input('');  % Valor de X teste "estimado"
          %---%disp('-----------------------------------------------------------------')
          %---%disp('                  Valor de Y para o X estimado:                ')
          %---%disp('-----------------------------------------------------------------')

          %---%printf('O valor de X estimado foi: %f\n',Xe)
          %---%printf('O valor de Y para o X estimado foi: %f\n',fcustom0(Xe)')
        %---%endif
        
        %---%disp('-----------------------------------------------------------------')
        %---%disp('             FIM DO AJUSTE NÃO LINEAR: CUSTOMIZADA               ')
        %---%disp('-----------------------------------------------------------------')
    
  case 5 %5)Gráfico de Dispersão
    disp('-----------------------------------------------------------------')
    disp('                     GRÁFICO DE DISPERSÃO                        ')
    disp('-----------------------------------------------------------------')
    disp('Insira os valores de X e Y para visualizar o diagrama de dispersão:')
    disp('INSIRA X, no formato de MATRIZ:')
    disp('EX Exponencial: [-1 -0.7 -0.4 -0.1 0.2 0.5 0.8 1]')
    disp('')
    
    x = input('');
    disp('-----------------------------------------------------------------')
    disp('INSIRA Y, no formato de MATRIZ:')
    disp('EX Exponencial: [36.547 17.264 8.155 3.852 1.820 0.860 0.406 0.246]')
    disp('')
    
    y = input('');
    
    scatter(x,y)
    
    disp('-----------------------------------------------------------------')
    disp('                ANALISE O GRÁFICO DE DISPERSÃO                   ')
    disp('-----------------------------------------------------------------')
    disp('Com base no gráfico de disperção, escolha qual função descreve   ')
    disp('melhor a dispersão dos valores inseridos.')
    disp('-----------------------------------------------------------------')
endswitch % FIM DO MENU DE AJUSTE LINEAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 12 %Plotar gráficos

%%%%%%%%%%%%%%%%%%%%%%% PLOTAR GRÁFICO %%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------------------------')
disp('                         PLOTAGEM DE GRÁFICO                       ')
disp('-------------------------------------------------------------------')
disp('Insira a função que deseja plotar:')
disp('EX: @(t) sin(t)')
disp('')
g = input('');

disp('Insira o intervalo que deseja visualizar do gráfico')
disp('EX: -4')
disp('Valor Inicial:')
disp('')
inicio = input('')
disp('EX: 6')
disp('Valor Final:')
disp('')
fim = input('')

%%FIM DO INPUT%%

disp('-------------------------------------------------------------------')
disp('                  GRÁFICO PLOTADO COM SUCESSO                      ')
disp('-------------------------------------------------------------------')
t = inicio:fim;

hold on
fplot(g,'-r')
hold off

grid on

legend('--Nenhuma Legenda aplicada--')
title('GRÁFICO PLOTADO COM SUCESSO')
xlabel('eixo x')
ylabel('eixo y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 13 %Informações adicionais sobre o código
disp('-------------------------------------------------------------------')
disp('                     INFORMAÇÕES ADICIONAIS                        ')
disp('-------------------------------------------------------------------')
disp('Utilize essa aba como forma de deixar possíveis erros que ocorreram')
disp('com você ou como um "Bloco de anotações" dentro do código. para que')
disp('você possa ter um acesso rápido a informações importantes.')
disp('')
disp('Tome cuidado na hora de aplicar funcões e equações dentro dos códigos.')
disp('Sempre utilize o @(variável da função) seguido de um espaço para ter')
disp('certeza que a equação seja inserida corretamente.')
disp('-------------------------------------------------------------------')
disp('            CRITÉRIOS DE PARADA QUE PODEM SER USADOS               ')
disp('-------------------------------------------------------------------')
disp('abs(f(p))<erro')
disp('abs(p - x0)<erro')
disp('abs(p - x0/p)<erro')
disp('-------------------------------------------------------------------')
disp('HUNGE KUTTA 3º e 4ºORDEM PROBLEMA!!')
disp('-------------------------------------------------------------------')
disp('Agradecimentos para: Mare, Galega, Borel, Caio...')
disp('Que ajudaram esse programa a se tornar possível, yeeeey Sz')
disp('Agradecimento especial para a professora Yara, que teve paciencia  ')
disp('com a gnt, o semestre de 2022.2 inteiro :)')
disp('-------------------------------------------------------------------')
endswitch %Sandwich

disp('Deseja usar o programa novamente?')
disp('1)SIM')
disp('2)NÃO')
disp('')
repe = input('');

disp('')

if repe == 1
  
  disp('Deseja limpar a janela de comandos antes de recomeçar?')
  disp('1)SIM')
  disp('2)NÃO')
  limpar = input('');
  
  if limpar == 1
    clear
    clc
    disp('                      CACHE APAGADO... T-T                         ')
    oneforall
  else
    disp('-------------------------------------------------------------------')
    disp('                      DADOS MANTIDOS... >.<                        ')
    oneforall
  endif
else
  disp('-------------------------------------------------------------------')
  disp('Até Mais, yeeey!!')
  disp('Qualquer coisa estou aqui...')
  disp('U.u   zZZZZzzzzzZzzzzZZz')
endif