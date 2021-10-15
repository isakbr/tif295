function alpha = chi_func(a,b,epsAu, eps0)
    if a == b
        chi = 2;
        alpha =4*pi*a^3*(epsAu-eps0)./(epsAu+chi*eps0);
    else
        xi = 1/sqrt(b^2/a^2-1);
        acos_term = ((xi*(xi^2+1))/2)*acos((xi^2-1)/(xi^2+1));
        chi = -1-2*1/(xi^2-acos_term);
        alpha = (epsAu-eps0)./(epsAu+chi*eps0)*(b^3*(1+chi)*xi^2+1)/(3*xi^2);
    end
end