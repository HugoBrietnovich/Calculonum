function[a,b,k]=menubissec(a,b,f,erro,modo)
%trocar função em outros casos%

%%Fazer função na janela como:%%
%% f = @(x) x .* log10 (x) - 1 %%
%% EX: [2,3]
%% Erro = 0.002

%%a e b são chutes iniciais%%

%%erro é dado%%

%usar @ usa f(x) e a entrada não tem aspas%
%se usar outro arquivo com a função, usa feval e o input da função vêm com aspas simples%

%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

disp('----------------------------------')

%disp('Defina qual método você deseja usar')
%disp('1) Bissecção')
%disp('2) Falsa posição')
%modo = input('');

disp('Defina os intervalos do "chute inicial"')

disp('Insira o valor de a:')
a = input('');

disp('Insira o valor de b:')
b = input('');

disp('Me informe a função que deseja aplicar o método da bissecção:')
disp('OBS: definir a função como: @(x) x .* log10 (x) - 1 POR EXEMPLO')
f = input('');

disp('Defina o erro a ser usado no critério de parada:')
erro = input('');

%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("  k        ak         bk          fak        fbk       pk       fpk       absfp")

if f(a)*f(b) > 0
   disp('"NAMORAL IRMÃO, TU TÁ NO TOP 3 PESSOAS MAIS BURRAS QUE EU CONHEÇO!!"')
   disp('Explica esse negócio pro tio... ~By Bill')
   error('Ta errado irmão, troca issae sobrinho ~By Bill!!')
else  
  k=0;
  
  while k<100 
      %if (modo = 1)
        p = (a+b)/2;
      %else (modo = 2)
        %p = a*abs(f(b))+b*abs(f(a))/abs(f(b))+abs(f(a));
      %endif 
    if abs(b-a)<erro
      disp('----------------------------------')
      printf('A raiz está entre %f e %f\n',a,b')
        %if (modo = 1)
          disp('Pelo método da Bissecção')
        %else modo = 2
          %disp('Pelo método da Falsa Posição')
        %endif
      break
      
    else  
      printf("%f  %f  %f  %f  %f  %f  %f  %f\n",k,a,b,f(a),f(b),p,f(p),abs(f(p)))
      if f(a)*f(p)<0
        b=p;
        
      else
        a=p;
        
      endif
      
      k=k+1;
    endif
    
  endwhile
endif
endfunction