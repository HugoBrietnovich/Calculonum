%   AJUSTE DE CURVAS   %
%  Vers�o 2.9kb  %

linear = menu('Ajuste de curvas:', '1)Ajuste Linear','2)Ajuste Polinomial', '3)Ajuste n�o linear','4)Ajuste Customizado', '5)Gr�fico de Dispers�o')
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

    %  SOMAT�RIOS  %
    sx = sum(x);
    sy = sum(y);
    sxQUAD = sum(x.^2);
    syQUAD = sum(y.^2);
    sxy = sum(x.*y);
    sx2 = (sum(x))^2;
    sy2 = (sum(y))^2;
    %              %

    disp('-----------------------------------------------------------------')
    disp('Deseja visualisar o valor dos som�t�rios?')
    disp('1)Sim')
    disp('2)Nao')
    disp('')
    versomats = input('');

    if versomats == 1
      disp('-----------------------------------------------------------------')
      printf('Somat�rio de x: %f\n',sx)
      printf('Somat�rio de y: %f\n',sy)
      printf('Somat�rio de x^2: %f\n',sxQUAD)
      printf('Somat�rio de x*y: %f\n',sxy)
    endif


    %  C�lculo dos Alfas, em somat�rio  %

    ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);

    ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);

    %  R^2  %
    disp('-----------------------------------------------------------------')
    disp('A equa��o da reta encontrada com os valores de Alfa 0 e Alfa 1 �:')
    disp('-----------------------------------------------------------------')
    printf('y = %f + (%f) * x\n',ALFA0, ALFA1)

    disp('')

    printf('Onde "Alfa 0" � igual a: %f\n',ALFA0) % Alfa 0 do exemplo: -20.078
    disp('e')
    printf('Onde "Alfa 1" � igual a: %f\n',ALFA1) % Alfa 1 do exemplo: 0.5276


    f0 = @(x) ALFA0 + ALFA1 * x;    % EQUA��O DA RETA COM PARES DADOS
    f1 = @(xe) ALFA0 + ALFA1 * xe;  % EQUA��O DA RETA COM X DE TESTE

    %  COEFICIENTE DE DETERMINA��O  %
    ymedio = sy/m;

    St = sum((y - ymedio).^2); % M�dia

    Sr = sum((y - f0(x)).^2);  % Curva ajustada

    RQUAD = (St-Sr)/St; % R^2 do exemplo: 0.8366

    RQuadComputacional = (((m*sxy)-(sx*sy))/(sqrt((m*sxQUAD)-sx2)*(sqrt((m*syQUAD)-sy2))))^2; % R^2 computacional funciona para o linear 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('-----------------------------------------------------------------')
    disp('                 Coeficiente de determina��o:                    ')
    disp('-----------------------------------------------------------------')

    printf('Valor do Y m�dio: %f\n',ymedio)
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
    disp('Deseja visualizar o diagrama de dispers�o?')
    disp('1)Sim')
    disp('2)N�o')
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
    disp('Insira o valor em X dos n�s (pares ordenados) em formato de MATRIZ')
    disp('EX: [300 350 400]')
    disp('')
    X = input('');

    disp('')
    
    disp('Insira o valor em Y dos n�s (pares ordenados) em formato de MATRIZ')
    disp('EX: [1.132 0.978 0.854]')
    disp('')
    Y = input('');
    
    disp('')
    
    disp('-----------------------------------------------------------------')
    disp('A ordem do polin�mio ser� calculada automaticamente, com base')
    disp('na quantidade de pares ordenados de n�s que foram inseridos.')
    n = (length(X)-1);
    
    disp('')
    
    printf('O polin�mio atual � de %f� ordem!\n',n)
    disp('-----------------------------------------------------------------')
    
    disp('Deseja alterar a ordem do polin�mio?')
    disp('1)Sim')
    disp('2)N�o')
    alterpoli = input('');
    
    if alterpoli == 1
      disp('Insira a nova ordem do polin�mio que ser� ajustado:')
      disp('')
      n = input('');
      if n < (length(X)-1)
        error('A ordem do polin�mio n�o pode ser menor que: %f',(length(X)-1))
      endif
    endif
    
    p = polyfit(X,Y,n);
    
    disp('-----------------------------------------------------------------')
    disp('O polin�mio obtido foi:')
    polyout(p)
    disp('-----------------------------------------------------------------')
    disp('                    FIM DO AJUSTE POLINOMIAL                     ')
    disp('-----------------------------------------------------------------')
    
  case 3 %  3)Ajuste n�o linear
    disp('Insira os valores de X e Y para visualizar o diagrama de dispers�o:')
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
    disp('                ANALISE O GR�FICO DE DISPERS�O                   ')
    disp('-----------------------------------------------------------------')
    disp('Com base no gr�fico de disper��o, escolha qual fun��o descreve   ')
    disp('melhor a dispers�o dos valores inseridos.')
    disp('-----------------------------------------------------------------')
    
    nlinear = menu('N�o lineares:', '1)Fun��o Exponencial', '2)Fun��o Logar�tmica', '3)Fun��o Potencial', '4)Fun��o Hiperb�lica')
    switch(nlinear)
      case 1 %  1)Fun��o Exponencial
     
        %  Exponencial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE N�O LINEAR: EXPONENCIAL                   ')
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
        
        disp('Deseja visualisar o valor dos som�t�rios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somat�rio de X: %f\n',sx)
          printf('Somat�rio de Y: %f\n',sy)
          printf('Somat�rio de X^2: %f\n',sxQUAD)
          printf('Somat�rio de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = exp(ALFA0);
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A fun��o exponencial ajustada para esses valores �:')
        printf('Y = %f * exp(%f * X)\n',a,b)
       
        disp('')

        printf('Onde "a" � igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" � igual a: %f\n',b) % b do exemplo:
       
        fexp0 = @(X) a * exp(b * X); % EQUA��O EXPONENCIAL COM PARES DADOS
        fexp1 = @(Xe) a * exp(b * Xe); % EQUA��O EXPONENCIAL COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINA��O  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % M�dia

        Sr = sum((Y - fexp0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gr�fico de Dispes�o                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        %disp('-----------------------------------------------------------------')
        %disp('                 Coeficiente de determina��o:                    ')
        %disp('-----------------------------------------------------------------')

        %printf('Valor do Y m�dio: %f\n',ymedio)
        %printf('Valor do St: %f\n',St)
        %printf('Valor do Sr: %f\n',Sr)

        %disp('-----------------------------------------------------------------')
        %disp('                      Valor do R qudradrado:                     ')
        %disp('-----------------------------------------------------------------')

        %printf('Valor do R^2: %f\n',RQUAD)
        
        disp('Deseja testar a fun��o encontrada?')
        disp('1)Sim')
        disp('2)N�o')
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
        disp('             FIM DO AJUSTE N�O LINEAR: EXPONENCIAL               ')
        disp('-----------------------------------------------------------------')
        
      case 2 %  2)Fun��o Logar�tmica %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE N�O LINEAR: LOGAR�TMICA                   ')
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
        
        disp('Deseja visualisar o valor dos som�t�rios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somat�rio de X: %f\n',sx)
          printf('Somat�rio de Y: %f\n',sy)
          printf('Somat�rio de X^2: %f\n',sxQUAD)
          printf('Somat�rio de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = ALFA1;
        b = exp(ALFA0/a);
        
        disp('-----------------------------------------------------------------')
        disp('A fun��o logar�tmica ajustada para esses valores �:')
        printf('Y = %f * log(%f * X)\n',a,b)
       
        disp('')

        printf('Onde "a" � igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" � igual a: %f\n',b) % b do exemplo:
       
        flogarit0 = @(X) a * log(b * X); % EQUA��O LOGAR�TMICA COM PARES DADOS
        flogarit1 = @(Xe) a * log(b * Xe); % EQUA��O LOGAR�TMICA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINA��O  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % M�dia

        Sr = sum((Y - flogarit0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gr�fico de Dispes�o                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a fun��o encontrada?')
        disp('1)Sim')
        disp('2)N�o')
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
        disp('             FIM DO AJUSTE N�O LINEAR: LOGAR�TMICA               ')
        disp('-----------------------------------------------------------------')
      
      case 3 %  3)Fun��o Potencial  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE N�O LINEAR: POTENCIAL                   ')
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
        
        disp('Deseja visualisar o valor dos som�t�rios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somat�rio de X: %f\n',sx)
          printf('Somat�rio de Y: %f\n',sy)
          printf('Somat�rio de X^2: %f\n',sxQUAD)
          printf('Somat�rio de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = exp(ALFA0);
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A fun��o potencial ajustada para esses valores �:')
        printf('Y = %f * (X^%f)\n',a,b)
       
        disp('')

        printf('Onde "a" � igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" � igual a: %f\n',b) % b do exemplo:
       
        fpoten0 = @(X) a * (X^b); % EQUA��O POTENCIAL COM PARES DADOS
        fpoten1 = @(Xe) a * (Xe^b); % EQUA��O POTENCIAL COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINA��O  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % M�dia

        Sr = sum((Y - fpoten0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gr�fico de Dispes�o                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a fun��o encontrada?')
        disp('1)Sim')
        disp('2)N�o')
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
        disp('              FIM DO AJUSTE N�O LINEAR: POTENCIAL                ')
        disp('-----------------------------------------------------------------')
      
      case 4 %  4)Fun��o Hiperb�lica  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE N�O LINEAR: HIPERB�LICA                   ')
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
        
        disp('Deseja visualisar o valor dos som�t�rios?')
        disp('1)Sim')
        disp('2)Nao')
        disp('')
        versomats = input('');

        if versomats == 1
          disp('-----------------------------------------------------------------')
          printf('Somat�rio de X: %f\n',sx)
          printf('Somat�rio de Y: %f\n',sy)
          printf('Somat�rio de X^2: %f\n',sxQUAD)
          printf('Somat�rio de X*Y: %f\n',sxy)
        endif
        
        ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        a = ALFA0;
        b = ALFA1;
        
        disp('-----------------------------------------------------------------')
        disp('A fun��o hiperb�lica ajustada para esses valores �:')
        printf('Y = %f + (%f/x)\n',a,b)
       
        disp('')

        printf('Onde "a" � igual a: %f\n',a) % a do exemplo:
        disp('e')
        printf('Onde "b" � igual a: %f\n',b) % b do exemplo:
       
        fhiperbo0 = @(X) a + (b/X); % EQUA��O HIPERB�LICA COM PARES DADOS
        fhiperbo1 = @(Xe) a + (b/Xe); % EQUA��O HIPERB�LICA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINA��O  %
        ymedio = sy/m;

        St = sum((Y - ymedio).^2); % M�dia

        Sr = sum((Y - fhiperbo0(X)).^2);  % Curva ajustada

        RQUAD = (St-Sr)/St; 
         
        disp('-----------------------------------------------------------------')
        disp('                   Novo Gr�fico de Dispes�o                      ')
        disp('-----------------------------------------------------------------')
        scatter(X,Y)
        
        disp('Deseja testar a fun��o encontrada?')
        disp('1)Sim')
        disp('2)N�o')
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
        disp('             FIM DO AJUSTE N�O LINEAR: HIPERB�LICA               ')
        disp('-----------------------------------------------------------------')
    endswitch
  
  case 4 %  4)Ajuste Customizado  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        disp('-----------------------------------------------------------------')
        disp('                AJUSTE N�O LINEAR: CUSTOMIZADO                   ')
        disp('-----------------------------------------------------------------')
        disp('')
        disp('OBSERVA��ES:')
        disp('OBS 1: Fa�a a an�lise da sua fun��o n�o linear e linearize-a.    ')
        disp('OBS 2: Identifique quem � o termo INDEPENDENTE e o DEPENDENTE.   ')
        disp('OBS 3: O termo DEPENDENTE estar� em fun��o do INDEPENDENTE.      ')
        disp('Por exemplo: Y(x) "Y em fun��o de X"                             ')
        disp('OBS 4: Idenfique quem ser� Y e X.                                ')
        disp('OBS 5: Identifique quem ser� (se existir) a e b.                 ')
        disp('OBS 6: EDITE esse trecho do c�digo, para sua lineariza��o customizada')
        disp('Basta pesquisar no editor: CUSTOMIZEAQUI')
        disp('')
        disp('-----------------------------------------------------------------')
        disp('                       EDITE O PROGRAMA                          ')
        disp('-----------------------------------------------------------------')
        
        % Tire do modo COMENT�RIO, as linhas com %---%
        
        %disp('Defina Y comparando com o y da equa��o que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('Y = log(y)')
        %disp('')
        
        %---%Y = %EDITE AQUI;
        
        %disp('Defina X comparando com o x da equa��o que vc quer linearizar:   ')
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
        
        %---%disp('Deseja visualisar o valor dos som�t�rios?')
        %---%disp('1)Sim')
        %---%disp('2)Nao')
        %---%disp('')
        %---%versomats = input('');

        %---%if versomats == 1
          %---%disp('-----------------------------------------------------------------')
          %---%printf('Somat�rio de X: %f\n',sx)
          %---%printf('Somat�rio de Y: %f\n',sy)
          %---%printf('Somat�rio de X^2: %f\n',sxQUAD)
          %---%printf('Somat�rio de X*Y: %f\n',sxy)
        %---%endif
        
        %---%ALFA0 = ((sxQUAD*sy)-(sxy*sx))/((m*sxQUAD)-sx2);
        %---%ALFA1 = ((m*sxy)-(sx*sy))/(m*sxQUAD-sx2);
        
        %---%disp('-----------------------------------------------------------------')
        %disp('Defina a comparando com o alfa0 da equa��o que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('a = exp(ALFA0)')
        %disp('')
        
        %---%a = %EDITE AQUI;
        
        %disp('Defina b comparando com o alfa1 da equa��o que vc quer linearizar:   ')
        %disp('EX: y = a*x^b "Potencial" --------- Y = alfa0 + alfa1*X          ')
        %disp('b = ALFA1')
        %disp('')
        
        %---%b = %EDITE AQUI;
        
        %---%disp('-----------------------------------------------------------------')
        %---%disp('A fun��o customuzada ajustada para esses valores �:')
        %---%printf('Y = EDITE AQUI %f e %f\n',a,b)
       
        %---%disp('')

        %---%printf('Onde "a" � igual a: %f\n',a) % a do exemplo:
        %---%disp('e')
        %---%printf('Onde "b" � igual a: %f\n',b) % b do exemplo:
       
        %fcustom0 = @(X) CUSTOMIZADA; % EQUA��O CUSTOMIZADA COM PARES DADOS
        %fcustom1 = @(Xe) CUSTOMIZADA; % EQUA��O CUSTOMIZADA COM X DE TESTE
        
        %  COEFICIENTE DE DETERMINA��O  %
        %---%ymedio = sy/m;

        %---%St = sum((Y - ymedio).^2); % M�dia

        %---%Sr = sum((Y - fcustom0(X)).^2);  % Curva ajustada

        %---%RQUAD = (St-Sr)/St; 
         
        %---%disp('-----------------------------------------------------------------')
        %---%disp('                   Novo Gr�fico de Dispes�o                      ')
        %---%disp('-----------------------------------------------------------------')
        %---%scatter(X,Y)
        
        %---%disp('Deseja testar a fun��o encontrada?')
        %---%disp('1)Sim')
        %---%disp('2)N�o')
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
        %---%disp('             FIM DO AJUSTE N�O LINEAR: CUSTOMIZADA               ')
        %---%disp('-----------------------------------------------------------------')
    
  case 5 %5)Gr�fico de Dispers�o
    disp('-----------------------------------------------------------------')
    disp('                     GR�FICO DE DISPERS�O                        ')
    disp('-----------------------------------------------------------------')
    disp('Insira os valores de X e Y para visualizar o diagrama de dispers�o:')
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
    disp('                ANALISE O GR�FICO DE DISPERS�O                   ')
    disp('-----------------------------------------------------------------')
    disp('Com base no gr�fico de disper��o, escolha qual fun��o descreve   ')
    disp('melhor a dispers�o dos valores inseridos.')
    disp('-----------------------------------------------------------------')
endswitch % FIM DO MENU DE AJUSTE LINEAR