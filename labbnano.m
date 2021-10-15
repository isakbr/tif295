%% Increasing discs
clf

lambda = linspace(100,900,700)
eps0 = 1
h = 4.135*1e-15 %eV/Hz
wpAu = 9 %eV
w = 3e8*2*pi;
c = 3e8
w_ev = h*c./(lambda*1e-9) %eV
epsAu = 1-(wpAu./w_ev).^2;
k = 2*pi./lambda
a = 20;
b = 20;
y_counter = 1
for i=2:2:21
    b = 20*(2+i/20)
    a = b
    xi = 1/sqrt(b^2/a^2-1);
    alpha = chi_func(a,b,epsAu,eps0)
    F = 1./(1-(2*j*k.^3.*alpha)/3-(k.^2.*alpha/b));
    Q = k.^4/(6*pi).*abs(alpha.*F).^2;
    Qmax = max(Q);
    %Q = Q/Qmax;
    y = ones(1,length(lambda));
    plot3(lambda,y*y_counter,Q,'color',[0 0.4470 0.7410],'Linewidth',2)
    hold on
    y_counter = y_counter + 1
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I (a.u.)$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;

%% Au
clf

wpAu = 7.4 %9.6 for au, 9.0 for ag % eV
epsAu = 1-(wpAu./w_ev).^2;
a = 20;
b = 20;
y_counter = 1
for i=2:2:21
    xi = 1/sqrt(b^2/a^2-1);
    alpha = chi_func(a,b,epsAu,eps0)
    F = 1./(1-(2*j*k.^3.*alpha)/3-(k.^2.*alpha/b));
    Q = k.^4/(6*pi).*abs(alpha.*F).^2;
    Qmax = max(Q);
    Q = Q/Qmax
    y = ones(1,length(lambda));
    plot3(lambda,y*y_counter,Q,'color',[0 0.4470 0.7410],'Linewidth',2)
    hold on
    y_counter = y_counter + 1
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('I (a.u)','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;