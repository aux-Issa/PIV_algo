% 対数則

%壁法則（τwは適当に代入）
kappa = 0.4;
B = 5.5;

% ucross = U_0 * sqrt(ro / tau);
ucross = uus(:,1)/ut * sqrt(ro / tauw);
% ycross = Yadjust / 1000 / nu * tauw ;
ycross = x(1:n,1,2) * sqrt(ro / tauw) / vis  ;
% ycross = ycross / 1000
% sqrt (tau / rho): 摩擦速度
% nu: 動粘度
% Yadjust：壁面ゼロからの座標

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
ucross2 = log(ycross1) / kappa + B

hold on
semilogx(ycross1,ucross1,'--');
semilogx(ycross1,ucross2,':');
xlim([0 1000]);
ylim([0 25]);
legend({'present','\slU^{+}= y^{+}','\slU^{+}=1/0.4lny^{+}+5.5'},'Location','southeast')

% todo: ，摩擦レイノルズ数で出してみる→対数速分布
% 渡辺の論文にuとv，Re応力の結果
% DNS（バルクレイノルズ数じゃなくて摩擦レイノルズ数）の結果　https://www.rs.tus.ac.jp/t2lab/db/
% 主流方向の速度分布を出してみる
% 壁面ませつけいすうを補正した値で壁面摩擦応力出してみる．