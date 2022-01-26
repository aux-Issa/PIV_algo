%壁法則（τwは適当に代入）
kappa = 0.4;
B = 5.5;

%ucross = U_0 * sqrt(rho / tau);
%ycross = Yadjust / 1000 / nu *sqrt (tau / rho);

figure(n)
%semilogx(ycross(1:162,1),ucross(1:162,1),'o');
grid on
grid minor
box on
xlabel('\sly^{+}','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\slU^{+}','FontName','Times','FontAngle','Italic','FontSize',20);
xlim([1 500])
ylim([0 25])
set(gca,'FontName','Times','FontSize',15)

ycross1 = logspace(-1,3);
ucross1 = ycross1;
ucross2 = log(ycross1) / kappa + B;
ucross3 = log(ycross1) / kappa1 + B1;

hold on
semilogx(ycross1,ucross1,'--');
semilogx(ycross1,ucross2,':');
semilogx(ycross1,ucross3,'-.');
hold off

legend({'present','\slU^{+}= y^{+}','\slU^{+}=1/0.4lny^{+}+5.5','\slU^{+}=1/0.21lny^{+}+0.6'},'Location','southeast')