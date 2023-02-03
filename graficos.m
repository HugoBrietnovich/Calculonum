function [f,g]=graficos(x)
f =@(x) x.^3 + sin(x)
g =@(x) x.^2 + cos(x)

hold on
fplot(f,'-or;x.^3 + sin(x);')

fplot(g,'-.+b;x.^2 + cos(x);')
hold off

grid on

legend('x.^3 + sin(x)','x.^2 + cos(x)')
title('seno e cosseno')
xlabel('eixo x')
ylabel('eixo y')
endfunction