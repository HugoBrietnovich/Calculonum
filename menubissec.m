function[a,b,k]=menubissec(a,b,f,erro,modo)
%trocar fun��o em outros casos%

%%Fazer fun��o na janela como:%%
%% f = @(x) x .* log10 (x) - 1 %%
%% EX: [2,3]
%% Erro = 0.002

%%a e b s�o chutes iniciais%%

%%erro � dado%%

%usar @ usa f(x) e a entrada n�o tem aspas%
%se usar outro arquivo com a fun��o, usa feval e o input da fun��o v�m com aspas simples%

%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

disp('----------------------------------')

%disp('Defina qual m�todo voc� deseja usar')
%disp('1) Bissec��o')
%disp('2) Falsa posi��o')
%modo = input('');

disp('Defina os intervalos do "chute inicial"')

disp('Insira o valor de a:')
a = input('');

disp('Insira o valor de b:')
b = input('');

disp('Me informe a fun��o que deseja aplicar o m�todo da bissec��o:')
disp('OBS: definir a fun��o como: @(x) x .* log10 (x) - 1 POR EXEMPLO')
f = input('');

disp('Defina o erro a ser usado no crit�rio de parada:')
erro = input('');

%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("  k        ak         bk          fak        fbk       pk       fpk       absfp")

if f(a)*f(b) > 0
   disp('"NAMORAL IRM�O, TU T� NO TOP 3 PESSOAS MAIS BURRAS QUE EU CONHE�O!!"')
   disp('Explica esse neg�cio pro tio... ~By Bill')
   error('Ta errado irm�o, troca issae sobrinho ~By Bill!!')
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
      printf('A raiz est� entre %f e %f\n',a,b')
        %if (modo = 1)
          disp('Pelo m�todo da Bissec��o')
        %else modo = 2
          %disp('Pelo m�todo da Falsa Posi��o')
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