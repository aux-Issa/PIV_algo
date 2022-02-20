%%
close all;  
clear;       
tic;




%%
%各数値データを入力
U_pito = 2.06;      %主流流速(m/s):ピトー管で計測
T = 20.0;           %実験温度(℃)
% L = 0.65;           %平板の長さ(m)
ls = 0.55;          %計測開始地点(m)
numfiles = 100;   %ファイル数
xn = 253;           %x方向のデータ数（列の数）
yn = 253;            %y方向のデータ数（行の数）
n = 1;

%%
%サザーランドの式による粘性係数の算出と動粘性係数の算出
T_0 = 0;             %基本温度(℃)
mu_0 = 1.72*10^-5;   %基本粘性係数(Pa*s)
C = 111;             %サザーランド係数(K)
P = 1013;            %気圧(hPa)
R = 2.87;            %乾燥空気の気体定数

T_0 = T_0 + 273.15;   %℃ -> K
T = T + 273.15;       %℃ -> K

mu = mu_0 * (T / T_0)^(3/2) * (T_0 + C) / (T + C);   %粘性係数(Pa*s)
rho = P / (R * T);                                    %空気の密度(kg/m^3)
nu = mu / rho;                                        %動粘性係数(m^2/s)


%%
%遷移の始まりの位置の数値計算予想
% Rex_min = 3.0*10^5;   %遷移レイノルズ数の最小値
% Rex_max = 5.0*10^5;   %遷移レイノルズ数の最大値
% 
% tr_s_min = Rex_min * nu / U_pito;       %遷移位置の最小値(m)
% tr_s_max = Rex_max * nu / U_pito;       %遷移位置の最大値(m)
% tr_s_avg = (tr_s_min + tr_s_max) / 2;   %平均値(m)


%%
%境界層厚さ（層流）
delta_laminar = 5 * sqrt(nu * ls / U_pito);
delta_turbulence = 0.37 * ls *(nu / U_pito / ls)^0.2;


%%
%ファイルからデータの読み取り
data = cell(numfiles,1);   %cell配列の大きさを定義

%     for j = 1:numfiles
%         if j-1 < 10
%            myfilename = sprintf('E:/doublepulselaser/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.00000%d.txt',j-1);
%         elseif j-1 < 100
%            myfilename = sprintf('E:/doublepulselaser/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.0000%d.txt',j-1);
%         elseif j-1 < 1000
%            myfilename = sprintf('E:/doublepulselaser/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.000%d.txt',j-1);
%         end
%     data{j,1} = readmatrix(myfilename);
%     end
    
    
    for j = 1:numfiles
        if j-1 < 10
           myfilename = sprintf('G:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.00000%d.txt',j-1);
        elseif j-1 < 100
           myfilename = sprintf('G:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.0000%d.txt',j-1);
        elseif j-1 < 1000
           myfilename = sprintf('G:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.000%d.txt',j-1);
        end
    data{j,1} = readmatrix(myfilename);
    end
% path要確認
%     for j = 1:numfiles
%         if j-1 < 10
%            myfilename = sprintf('F:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.00000%d.txt',j-1);
%         elseif j-1 < 100
%            myfilename = sprintf('F:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.0000%d.txt',j-1);
%         elseif j-1 < 1000
%            myfilename = sprintf('F:/APG/u2.0ms_range50_100images_validation/u2.0ms.6s3cw2c8.000%d.txt',j-1);
%         end
%     data{j,1} = readmatrix(myfilename);
%     end


%%
%x,y座標データの読み取り
    X_swap = zeros(yn,xn);   %ゼロ行列として定義
    Y_swap = zeros(yn,xn);   %ゼロ行列として定義

    x = data{1,1}(:,5);           %x座標
    X = reshape(x,yn,xn);         %x座標をyn*xn行列に変換
    X = X.';
    X = X + 525;
%      for i = 1:yn
%          X_swap(:,i) = X(:,yn+1-i);
%      end
%      X(:,:) = X_swap(:,:);
%     X = sort(X,2);          %x座標を左右反転（左原点）

    y = data{1,1}(:,6);                  %y座標
    Y = reshape(y,yn,xn);                %y座標をyn*xn行列に変換
    Y = Y.';
    Y = sort(Y,'descend');               %y座標を降順にする
    Yadjust = Y - 19;
%     for j = 1:xn
%         Y_swap(:,j) = Y(:,xn+1-j);
%     end
%     Y(:,:) = Y_swap(:,:);                %y座標を左右反転（左原点）（この操作はやらなくても変わらないはず）


%%
%速度u,vの読み取り
u = cell(numfiles,1);   %cell配列の大きさを定義
U = cell(numfiles,1);   %cell配列の大きさを定義
v = cell(numfiles,1);   %cell配列の大きさを定義
V = cell(numfiles,1);   %cell配列の大きさを定義
U_swap = zeros(yn,xn);     %ゼロ行列として定義
V_swap = zeros(yn,xn);     %ゼロ行列として定義

    for j = 1:numfiles
        u{j,1}(:) = data{j,1}(:,9);               %x方向の速度成分u
        u{j,1}(:) = -1 * u{j,1}(:);               %uを正の値に変換
        U{j,1}(:,:) = reshape(u{j,1}(:),yn,xn);   %uをyn*xn行列に変換
        U{j,1}(:,:) = U{j,1}(:,:).';
 
        v{j,1}(:) = data{j,1}(:,10);               %x方向の速度成分u
        V{j,1}(:,:) = reshape(v{j,1}(:),yn,xn);   %uをyn*xn行列に変換
        V{j,1}(:,:) = V{j,1}(:,:).';        
        
        for i = 1:yn
            U_swap(i,:) = U{j,1}(yn+1-i,:);
            V_swap(i,:) = V{j,1}(yn+1-i,:);
        end
        U{j,1}(:,:) = U_swap(:,:);
        V{j,1}(:,:) = V_swap(:,:);
        
        for i = 1:xn
            U_swap(:,i) = U{j,1}(:,xn+1-i);
            V_swap(:,i) = V{j,1}(:,xn+1-i);
        end
        U{j,1}(:,:) = U_swap(:,:);
        V{j,1}(:,:) = V_swap(:,:);
        
    end


%%
%平均主流流速を算出
U_sum = zeros(yn,xn);
U_tavg_sum = zeros(yn,1);
V_sum = zeros(yn,xn);
V_tavg_sum = zeros(yn,1);

    for j = 1:numfiles
        U_sum = U_sum + U{j,1}(:,:);
        V_sum = V_sum + V{j,1}(:,:);
    end
    U_tavg = U_sum / numfiles;
    V_tavg = V_sum / numfiles;
    U_tavg1 = round(U_tavg,2);
    V_tavg1 = round(V_tavg,2);
    
    for j = 33:63
        U_tavg_sum(:,1) = U_tavg_sum(:,1) + U_tavg1(:,j);
        V_tavg_sum(:,1) = V_tavg_sum(:,1) + V_tavg1(:,j);
    end
    for j = 76:118
        U_tavg_sum(:,1) = U_tavg_sum(:,1) + U_tavg1(:,j);
        V_tavg_sum(:,1) = V_tavg_sum(:,1) + V_tavg1(:,j);
    end
        U_0 = U_tavg_sum / 74;
        V_0 = V_tavg_sum / 74;

figure(n)
contourf(X(1:162,:),Yadjust(1:162,:),U_tavg(1:162,:),12)
xlabel('\slx \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylim([0 30]);
c = colorbar;
c.Label.String = '\slU_{0} \rm(m/s)';
c.Label.FontSize = 20;
c.Label.FontName = 'Times';
caxis([0 2.2]);
n = n + 1;

%瞬時用
% for i = 91:100
%     figure(n)
%     contourf(X(1:162,:),Yadjust(1:162,:),U{i,1}(1:162,:),7)
%     xlabel('\slx \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
%     ylabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
%     ylim([0 30]);
%     c = colorbar;
%     c.Label.String = '\slU_{0} \rm(m/s)';
%     c.Label.FontSize = 20;
%     c.Label.FontName = 'Times';
%     caxis([0 2.2]);
%     n = n + 1;
% end

figure(n)
sz = 20;
scatter(U_0(1:162,1),Yadjust(1:162,1),sz);
grid on
grid minor
box on
xlabel('\slu \rm(m/s)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
set(gca,'FontName','Times','FontSize',15)
n = n + 1;

% figure(n)
% sz = 20;
% scatter(V_0,Y(:,1),sz);
% grid on
% grid minor
% box on
% xlabel('\slu / \slU_{0}','FontName','Times','FontAngle','Italic','FontSize',20);
% ylabel('\sly','FontName','Times','FontAngle','Italic','FontSize',20);
% set(gca,'FontName','Times','FontSize',15)
% n = n + 1;


%壁法則（τwは適当に代入）
tau = 0.009;
kappa = 0.4;
kappa1 = 0.21;
B = 5.5;
B1 = 0.6;

ucross = U_0 * sqrt(rho / tau);
ycross = Yadjust / 1000 / nu *sqrt (tau / rho);

figure(n)
semilogx(ycross(1:162,1),ucross(1:162,1),'o');
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
n = n + 1;


%%
%u'RMS,v'RMS,レイノルズせん断応力
u_henndou = cell(numfiles,1);    %cell配列の大きさを定義
u_henndouS = cell(numfiles,1);   %cell配列の大きさを定義
u_henndouS_sum = zeros(yn,xn);
v_henndou = cell(numfiles,1);    %cell配列の大きさを定義
v_henndouS = cell(numfiles,1);   %cell配列の大きさを定義
v_henndouS_sum = zeros(yn,xn);
taure1 = cell(numfiles,1);
taure_sum = zeros(yn,xn);

    for j = 1:numfiles
        u_henndou{j,1}(:,:) = U{j,1}(:,:) - U_tavg(:,:);
        u_henndouS{j,1} = u_henndou{j,1} .* u_henndou{j,1};
        u_henndouS_sum = u_henndouS_sum + u_henndouS{j,1}(:,:);
        
        v_henndou{j,1}(:,:) = V{j,1}(:,:) - V_tavg(:,:);
        v_henndouS{j,1} = v_henndou{j,1} .* v_henndou{j,1};
        v_henndouS_sum = v_henndouS_sum + v_henndouS{j,1}(:,:);
        
        taure1{j,1}(:,:) = u_henndou{j,1}(:,:) .* v_henndou{j,1}(:,:);
        taure_sum = taure_sum + taure1{j,1}(:,:);
    end
    u_henndouMS = u_henndouS_sum / numfiles;
    u_henndouRMS = sqrt(u_henndouMS);

    v_henndouMS = v_henndouS_sum / numfiles;
    v_henndouRMS = sqrt(v_henndouMS);
    
    taure_tavg = taure_sum / numfiles;
    
    urms_sum = zeros(yn,1);
    vrms_sum = zeros(yn,1);
    taure_tavg_sum = zeros(yn,1);
    for i = 33:63
        urms_sum = urms_sum + u_henndouRMS(:,i);
        vrms_sum = vrms_sum + v_henndouRMS(:,i);
        taure_tavg_sum = taure_tavg_sum + taure_tavg(:,i);
    end
    for i = 76:118
        urms_sum = urms_sum + u_henndouRMS(:,i);
        vrms_sum = vrms_sum + v_henndouRMS(:,i);
        taure_tavg_sum = taure_tavg_sum + taure_tavg(:,i);
    end
    urms = urms_sum / 74;
    vrms = vrms_sum / 74;
    taure2 = taure_tavg_sum / 74;
    taure = -1 * rho * taure2;
    
    uvrms = u_henndouRMS .* v_henndouRMS;
    Ruv = -1 * taure_tavg ./ uvrms;
    
    Ruv_sum = zeros(yn,1);
%指定した範囲or全体での流れ方向平均
%     for i = 33:63
%         Ruv_sum = Ruv_sum + Ruv(:,i);
%     end
%     for i = 76:118
%         Ruv_sum = Ruv_sum + Ruv(:,i);
%     end
%     Ruv_xavg = Ruv_sum / 74;
    
    for i = 1:xn
        Ruv_sum = Ruv_sum + Ruv(:,i);
    end
    Ruv_xavg = Ruv_sum / xn;

figure(n)
sz = 20;
scatter(Yadjust(1:162,1),urms(1:162,1),sz);
grid on
grid minor
box on
xlabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\slu_{rms}^{,} \rm(m/s)','FontName','Times','FontAngle','Italic','FontSize',20);
set(gca,'FontName','Times','FontSize',15)
n = n + 1;

figure(n)
sz = 20;
scatter(Yadjust(1:162,1),vrms(1:162,1),sz);
grid on
grid minor
box on
xlabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\slv_{rms}^{,} \rm(m/s)','FontName','Times','FontAngle','Italic','FontSize',20);
set(gca,'FontName','Times','FontSize',15)
n = n + 1;

figure(n)
sz = 20;
scatter(Yadjust(1:162,1),taure(1:162,1),sz);
grid on
grid minor
box on
str = '$$ \tau = - \rho \overline{u^{,}v^{,}} $$';
xlabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel(str,'Interpreter','latex','FontSize',20);
set(gca,'FontName','Times','FontSize',15)
n = n + 1;

figure(n)
contourf(X(1:162,:),Yadjust(1:162,:),Ruv(1:162,:),6)
xlabel('\slx \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
 ylim([0 35]);
c = colorbar;
c.Label.String = '\slR_{u^{,}v^{,}}';
c.Label.FontSize = 20;
c.Label.FontName = 'Times';
caxis([-1 1]);
n = n + 1;

figure(n)
sz = 20;
scatter(Ruv_xavg(1:162,1),Yadjust(1:162,1),sz);
grid on
grid minor
box on
xlabel('\slR_{u^{,}v^{,}}','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\sly \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
set(gca,'FontName','Times','FontSize',15)
n = n + 1;


%%
%Cp分布
cp_ysum = zeros(1,xn);

for i = 16:36
    cp_ysum = U_tavg(i,:) + cp_ysum;
end
cp_yavg = cp_ysum / 21;
% cp_yavg = round(cp_yavg,2);

cp1 = cp_yavg / cp_yavg(1,76);
cp2 = cp1.^2;
cp = 1 - cp2;

% figure(n)
% sz = 20;
% scatter(X(1,1:194),cp(1,1:194),sz);
% grid on
% grid minor
% box on
% xlabel('\slx \rm(mm)','FontName','Times','FontAngle','Italic','FontSize',20);
% ylabel('\slC_{p}','FontName','Times','FontAngle','Italic','FontSize',20);
% set(gca,'FontName','Times','FontSize',15)
% ylim([-0.2 0.2]);
% n = n + 1;


%%
%ブラジウス解との比較
% Parameters of Blasius Equation
U_inf = 1;
L = 10;
A = sqrt(nu/U_inf);
h = 0.01;

% Numerical Solution of Blasius Equation Using Runge-Kutta
f1 = @(x, y1, y2, y3) y2;
f2 = @(x, y1, y2, y3) y3;
f3 = @(x, y1, y2, y3) -y1*y3;
eta = 0:h:10;
x = 0:h:10;
y1(1) = 0;
y2(1) = 0;
y3(1) = 0.4696;
for i = 1:(length(eta)-1)
a = h.*[f1(eta(i), y1(i), y2(i), y3(i)), f2(eta(i), y1(i), y2(i), y3(i)), f3(eta(i), y1(i), y2(i), y3(i))];
b = h.*[f1(eta(i), y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f2(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f3(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2)];
c = h.*[f1(eta(i), y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f2(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f3(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2)];
d = h.*[f1(eta(i), y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f2(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f3(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3))];
y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end

% Plotting and Visualization
% figure(1)
% plot(y1, eta, y2, eta, y3, eta, 'LineWidth', 2)
% xlim([0 2])
% title('Solution of Blasius eqution', 'FontSize', 14);
% xlabel('f, f'' and f''''', 'FontSize', 20);
% ylabel('\eta', 'FontSize', 20);
% grid on
% Legend1 = {'f(\eta)', 'f''(\eta)', 'f''''(\eta)'};
% legend(Legend1, 'FontSize', 14);

eta1 = sqrt(U_pito / nu / 0.55);
eta2 = Yadjust(:,1) * eta1 / 1000;
ublasius = U_0 / U_pito;

eta4sum = zeros(yn,1);
for i = 33:63
    eta3 = sqrt(U_pito / nu /X(1,i));
    eta4 = Yadjust(:,1) * eta3 / 1000;
    eta4sum = eta4sum + eta4;
end
for i = 76:118
    eta3 = sqrt(U_pito / nu /X(1,i));
    eta4 = Yadjust(:,1) * eta3 / 1000;
    eta4sum = eta4sum + eta4;
end
eta5 = eta4sum / 74;

% b = Y(:,1) ./ delta_turbulence;
% a = U_pito .* (b .^ (1/7));

figure(n)
sz = 20;
scatter(ublasius(1:162,1),eta2(1:162,1),sz);  %eta2は時間平均したもの　eta5はそれぞれのxで計算し平均したもの
grid on
grid minor
box on
xlabel('\slu / \slU_{0}','FontName','Times','FontAngle','Italic','FontSize',20);
ylabel('\eta = \sly(\slU_{0}/\nu\slx)^{1/2}','FontName','Times','FontAngle','Italic','FontSize',20);
xlim([0 1.0]);
set(gca,'FontName','Times','FontSize',15)

hold on
plot(y2,eta);
hold off

legend({'present','Blasius boundary layer'},'Location','northwest')

n = n + 1;


%%
%排除厚さと運動量厚さの計算，最後に形状係数の算出
Yadjust_T = Yadjust(1:162,1)';   %境界層厚さを求めるときの積分範囲であるy座標
Yadjust_T(1,163) = 0;            %壁面上の座標を追加:y=0

UH = U_0(1:162,1);
UH(163,1) = 0;

Umu = UH - UH(1,1);
utUmu = UH .* Umu;
Umu_T = Umu';
utUmu_T = utUmu';

blt_dis_T = trapz(Yadjust_T,Umu_T,2);     %(U-u)を台形積分
blt_mom_T = trapz(Yadjust_T,utUmu_T,2);   %u(U-u)を台形積分

blt_dis_TT = blt_dis_T';
blt_dis = blt_dis_TT / UH(1,1);    %排除厚さのx方向分布

UH2 = UH(1,1)^2;
blt_mom_TT = blt_mom_T';
blt_mom = blt_mom_TT / UH2;   %運動量厚さのx方向分布

H = blt_dis / blt_mom;            %形状係数の算出（x方向分布）

% %排除厚さδ*分布図を出力
% figure(n);                            %n枚目の図
% sz = 20;                              %プロットサイズ変更
% scatter(X(1,:),blt_dis{setn,1},sz);   %散布図
% xlabel('\slx\rm/L');                  %x軸ラベル
% ylabel('\delta^*');                   %y軸ラベル
% title('Displacement thickness');      %図の名前
% n = n + 1;                            %図の番号調整
% 
% %排除厚さθ分布図を出力
% figure(n);                            %n枚目の図
% sz = 20;                              %プロットサイズ変更
% scatter(X(1,:),blt_mom{setn,1},sz);   %散布図
% xlabel('\slx\rm/L');                  %x軸ラベル
% ylabel('\theta');                      %y軸ラベル
% title('Momentum thickness');          %図の名前
% n = n + 1;                            %図の番号調整
% 
%形状係数H分布図を出力
% figure(n);                      %n枚目の図
% sz = 20;                        %プロットサイズ変更
% scatter(X(1,:),H{setn,1},sz);   %散布図
% xlabel('\slx\rm/L');            %x軸ラベル
% ylabel('\slH');                 %y軸ラベル
% title('Shape factor');          %図の名前
% n = n + 1;                      %図の番号調整


%%
toc;   %計算時間計測終了

