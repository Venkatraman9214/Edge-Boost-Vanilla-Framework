function ODE=ODE_ICME(INI,lambda)
beta = lambda;
psi = 0.03;
[t,ODE]=ode45('ICME',[0,300],INI,[],beta,psi);
end
