function vals = matsolv(pars, x, Delta)

A1=pars(1);
b1=pars(2);
LW1=pars(3);
A2=pars(4);
b2=pars(5);
LW2=pars(6);



off=x;



p=A1*1./(1+((off-b1)./(0.5.*LW1)).^2)+ A2*1./(1+((off-b2)./(0.5.*LW2)).^2);

 vals=p;
end
