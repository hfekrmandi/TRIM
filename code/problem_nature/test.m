function test
close all
rho_hat  = 100;
for c=1:540
    p(c) = erlang_(c,rho_hat);
end
plot(p)


end
function p= erlang_(c,rho_hat)
rho = rho_hat;
den = 0;
for i=1:c
    den = den + (rho^i)/factorial(i);
end
num = (rho^c)/factorial(c);
p = num/den;
end