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