function value = myOwenT(h,a)
% this is self-defined Owen's T funciton, which is used when fitting a
% skewed gaussian distribution
fun=@(x) exp(-h.^2.*(1+x.^2)/2)./(1+x.^2);
F=@(y)integral(fun,0,y,'ArrayValued',true);
value = F(a)/2/pi;
end