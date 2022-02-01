% 対数則

%壁法則（τwは適当に代入）
kappa = 0.4;
B = 5.5;
fixed_tauw = 0.374;
u_tau = sqrt(fixed_tauw / ro)  %　摩擦速度
ucross = ums(:,1) / u_tau;
% 壁からの距離
y_coord = x(1:n,1,2);
ycross = y_coord * u_tau / vis;
semilogx(ycross,ucross,'o');
% ycross = ycross / 1000
% sqrt (tau / rho): 摩擦速度
% nu: 動粘度
% Yadjust：壁面ゼロからの座標

figure(n)
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
semilogx(ycross1,ucross1,'--');
hold on
semilogx(ycross,ucross,'o');
semilogx(ycross1,ucross2,':');
xlim([0 1000]);
ylim([0 25]);
legend({'\slU^{+}=1/0.4lny^{+}+5.5','present','\slU^{+}= y^{+}'},'Location','southeast')
hold off

% todo: ，摩擦レイノルズ数で出してみる→対数速分布
% 渡辺の論文にuとv，Re応力の結果
% DNS（バルクレイノルズ数じゃなくて摩擦レイノルズ数）の結果　https://www.rs.tus.ac.jp/t2lab/db/
% 主流方向の速度分布を出してみる
% 壁面ませつけいすうを補正した値で壁面摩擦応力出してみる．
% ums(:,1): uの時空間平均
