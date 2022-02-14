%%
close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
tic;                                % start of measuring the elapsed time


%%
% parameter part
numfiles = 500;                     % total number of file                      
m        = 127;                      % number of x direction
n        = 108;                      % number of y direction
ub       = 0.335;                    % bulk velocity [m/s]
% 補正前の壁面せん断応力                  
% tauw     = 0.23;

% excelから算出した補正後の壁面せん断応力
tauw     = 0.095331852;                   % wall shear stress [N/m^2]
ro       = 997.05;                  % density [kg/m^3]
ut       = sqrt(tauw/ro);           % friction velocity [m/s] 
vis      = 0.000000893;             % kinetic viscosity [m^2/s]
h        = 0.02;                    % chanel half height [m]
um       = zeros(n,m,2);            % matlix for velocity 
uu       = zeros(n,m,3);            % matlix for velocity fluctuation

%  input vec. 
%                  wall
%                   m
%              |---------|ch
%             n|  Flow�� |
%           ��y|---------|
%             o �� x  
% It is noted that the normal direction is reversed in the calculation owting to the wall location.
%%
% reading part_coordinate_data
loadname   = ('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/water_experiment.6vari2x3/coordinate_32px_v2.6vas9dlb.000000.dat');                                        % inputting file name
coord      = importdata(loadname);                                 % making structure of above
x(:,:,1)   = rot90(reshape(coord.data(:,1),[m n]));                % x data [mm]
x(:,:,2)   = rot90(reshape(coord.data(:,2),[m n]));                % y data [mm]
x(:,:,2)   = flip(x(:,:,2),1);                                     % y data [mm] with flip
x          = x/1000;                                               % xy data [mm]��[m]


%%
% reading part_velocity_data_average
for k =0:numfiles-1
	if k<10                                                        % changing file name to read data
	 file_name = sprintf('velocity_32px.6vayg6cd.00000%d.dat',k);
	elseif k<100
	 file_name = sprintf('velocity_32px.6vayg6cd.0000%d.dat',k);
	elseif k<1000
	 file_name = sprintf('velocity_32px.6vayg6cd.000%d.dat',k);
	elseif k<10000 
	 file_name = sprintf('velocity_32px.6vayg6cd.00%d.dat',k);
    end
    myfilename = sprintf('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/80ppm_experiment.6vaxq7zy/velocity/%s',file_name);
    mydata     = importdata(myfilename);
    u(:,:,1)   = rot90(reshape(mydata.data(:,1),[m n]));            % u
    u(:,:,2)   = -rot90(reshape(mydata.data(:,2),[m n]));           % v

    um(:,:,1)  = um(:,:,1)+u(:,:,1)/double(numfiles);               % time ave. u
    um(:,:,2)  = um(:,:,2)+u(:,:,2)/double(numfiles);               % time ave. v
end

    ums(:,1)  = mean(um(:,:,1),2);                                  % emsemble ave. u
    ums(:,2)  = mean(um(:,:,2),2);                                  % emsemble ave. u

% %%
% % 2D-figure
figure;
newpoints = m*3;
[xq,yq] = meshgrid(...                                              % making fine mesh
               linspace(x(1,1,1),x(1,m,1),newpoints ),...
               linspace(x(1,1,2),x(n,1,2),newpoints )...
             );
Umxy = interp2(x(:,:,1),x(:,:,2),um(:,:,1),xq,yq,'spline');         % interpolation
[cc,hc]=contourf(xq/h,yq/h,Umxy/ub,8);                             % making isoline
%[cc,hc]=contourf(xq*ut/vis,yq*ut/vis,Umxy/ut,12);                   % making isoline
hc.TextStep = 0.4;                                                  % interval isoline
hc.ShowText = 'off';                                                % isoline text off
colormap('jet');                                                    % color type of 'jet'
box on;                                                             % making flame of figure
% xlim([0 4]);                                                       % range of x
% ylim([0 0.5]);                                                     % range of y
set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
ax = gca;
ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
xlabel('{\it x}/{\it h}','FontSize',20,'Interpreter','latex');      % xlabel, its size, and type
ylabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');      % ylabel, its size, and type
%xlabel('{\it x}^+','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
%ylabel('{\it y}^+','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
c = colorbar;                                                       % making color bar
c.Limits = [0 1.2];                                                 % range of colorbar
c.FontSize = 18;                                                    % font size of scale of color bar
c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
c.Label.Interpreter = 'latex';                                      % font type of label of color bar
c.Label.String = '$${\it \overline{u}}/{\it U_b}$$';                % label of color bar
%c.Label.String = '$${\it \overline{u^+}}$$';                       % label of color bar
c.Label.FontSize = 20;                                              % font size of label of color bar
saveas(gcf,'Umxy','png');                                           % output
% 
% %%
% 1D-figure
figure;
% p = plot(decimate(x(1:n,1,2),2)/h,decimate(ums(:,1),2)/ub,'ko');
p = plot(x(1:n,1,2)/h,ums(:,1)/ub,'ko');
p.LineWidth = 1.2;
p.MarkerSize = 10;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.1 1 1];
box on;
%axis equal tight;                                                 
xticks([0:0.2:1]);
yticks([0:0.2:1.4]);
xlim([0 1]);
ylim([0 1.3]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
ylabel('$${\it \overline{u}}/{\it U_b}$$','FontSize',20,'Interpreter','latex');
saveas(gcf,'Um-y','png');

%%
% reading part_velocity_data_fluctuation
for k =0:numfiles-1
	if k<10                                                         % changing file name to read data
	 file_name = sprintf('velocity_32px.6vayg6cd.00000%d.dat',k);
	elseif k<100
	 file_name = sprintf('velocity_32px.6vayg6cd.0000%d.dat',k);
	elseif k<1000
	 file_name = sprintf('velocity_32px.6vayg6cd.000%d.dat',k);
	elseif k<10000 
	 file_name = sprintf('velocity_32px.6vayg6cd.00%d.dat',k);
    end
    myfilename = sprintf('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/80ppm_experiment.6vaxq7zy/velocity/%s',file_name);
    mydata     = importdata(myfilename);
    u(:,:,1)   = rot90(reshape(mydata.data(:,1),[m n]));            % u
    u(:,:,2)   = -rot90(reshape(mydata.data(:,2),[m n]));           % v
    % 変動を算出
    uf(:,:,1)  = u(:,:,1)-ums(:,1);
    uf(:,:,2)  = u(:,:,2)-ums(:,2);
    uf(:,:,3)  = -uf(:,:,1).*uf(:,:,2);
    % 変動の二乗平均
    uu(:,:,1)  = uu(:,:,1)+uf(:,:,1).*uf(:,:,1)/double(numfiles);
    uu(:,:,2)  = uu(:,:,2)+uf(:,:,2).*uf(:,:,2)/double(numfiles);
    uu(:,:,3)  = uu(:,:,3)+uf(:,:,3)/double(numfiles);
    % 変動の三乗平均(skewness)
    base = zeros(n,m,3);
    uu_Skewness(:,:,1)  = base(:,:,1)+uf(:,:,1).*uf(:,:,1).*uf(:,:,1)/double(numfiles);
    uu_Skewness(:,:,2)  = base(:,:,2)+uf(:,:,2).*uf(:,:,2).*uf(:,:,2)/double(numfiles);

    % 変動の四乗平均(flatness)
    uu_Flatness(:,:,1)  = base(:,:,1)+uf(:,:,1).*uf(:,:,1).*uf(:,:,1).*uf(:,:,2)/double(numfiles);
    uu_Flatness(:,:,2)  = base(:,:,2)+uf(:,:,2).*uf(:,:,2).*uf(:,:,2).*uf(:,:,2)/double(numfiles);
    % uu_Skewness(:,:,3)  = uu(:,:,3)+uf(:,:,3)/double(numfiles);
    U_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,1),xq,yq,'spline');         % interpolation
    V_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,2),xq,yq,'spline');         % interpolation
    U_dash_times_V_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,3),xq,yq,'spline');         % interpolation

% uの変動を分布図に

    % if k < 15
    % if k == 6 || k == 8|| k == 14
    %     [cc_solution_dash,hc_solution_dash]=contourf(xq*ut/vis,yq*ut/vis,U_dash/ut,8);                   % making isoline
    %     hc_solution_dash.TextStep = 0.4;                                                  % interval isoline
    %     hc_solution_dash.ShowText = 'off';                                                % isoline text off
    %     colormap('jet');                                                    % color type of 'jet'
    %     box on;                                                             % making flame of figure
    %     xlim([0 150]);                                                       % range of x
    %     ylim([0 60]);                                                     % range of y
    %     set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
    %     ax = gca;
    %     ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
    %     xlabel('${\it x^+}$','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
    %     ylabel('${\it y^+}$','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
    %     c = colorbar;                                                       % making color bar
    %     % c.Limits = [-0.4 0.4];                                                 % range of colorbar
    %     c.Limits = [-6 4];                                                 % range of colorbar
    %     c.FontSize = 18;                                                    % font size of scale of color bar
    %     c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
    %     c.Label.Interpreter = 'latex';                                      % font type of label of color bar
    %     % c.Label.String = '$${\it u^\prime}/{\it U_b}$$';                % label of color bar
    %     c.Label.String = '$${\it u^\prime}^+$$';                             % label of color bar
    %     c.Label.FontSize = 20;                                              % font size of label of color bar
    %     fig_name = sprintf('U_dash_solution%d',k);
    %     saveas(gcf,fig_name ,'png'); 
        
    % end 
    % % if k < 15
    % if k == 6 || k == 8|| k == 14
    %     [cc_solution_dash,hc_solution_dash]=contourf(xq*ut/vis,yq*ut/vis,V_dash/ut,8);                   % making isoline
    %     hc_solution_dash.TextStep = 0.4;                                                  % interval isoline
    %     hc_solution_dash.ShowText = 'off';                                                % isoline text off
    %     colormap('jet');                                                    % color type of 'jet'
    %     box on;                                                             % making flame of figure
    %     xlim([0 150]);                                                       % range of x
    %     ylim([0 60]);                                                     % range of y
    %     set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
    %     ax = gca;
    %     ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
    %     xlabel('${\it x^+}$','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
    %     ylabel('${\it y^+}$','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
    %     c = colorbar;                                                       % making color bar
    %     % c.Limits = [-0.4 0.4];                                                 % range of colorbar
    %     c.Limits = [-3 3];
    %     c.FontSize = 18;                                                    % font size of scale of color bar
    %     c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
    %     c.Label.Interpreter = 'latex';                                      % font type of label of color bar
    %     % c.Label.String = '$${\it v^\prime}/{\it U_b}$$';                % label of color bar
    %     c.Label.String = '$${\it u^\prime}^+$$';                             % label of color bar
    %     c.Label.FontSize = 20;                                              % font size of label of color bar
    %     fig_name = sprintf('V_dash_solution%d',k);
    %     saveas(gcf,fig_name ,'png'); 
        
    % end  
    % レイノルズ剪断応力の瞬時場
    % if k < 15
    % if k == 6 || k == 8|| k == 14
    %     [cc_solution_ReStress,hc_solution_ReStress]=contourf(xq*ut/vis,yq*ut/vis,U_dash_times_V_dash/(ut*ut),16);                   % making isoline
    %     hc_solution_ReStress.TextStep = 0.4;                                                  % interval isoline
    %     hc_solution_ReStress.ShowText = 'off';                                                % isoline text off
    %     colormap('jet');                                                    % color type of 'jet'
    %     box on;                                                             % making flame of figure
    %     xlim([0 150]);                                                       % range of x
    %     ylim([0 60]);                                                     % range of y
    %     set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
    %     ax = gca;
    %     ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
    %     xlabel('${\it x^+}$','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
    %     ylabel('${\it y^+}$','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
    %     c = colorbar;                                                       % making color bar
    %     % c.Limits = [-0.05 0.05];                                                 % range of colorbar
    %     c.Limits = [-5 15];
    %     c.FontSize = 18;                                                    % font size of scale of color bar
    %     c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
    %     c.Label.Interpreter = 'latex';                                      % font type of label of color bar
    %     c.Label.String = '$$-\overline{{\it u^\prime}^+{\it v^\prime}^+}$$';                % label of color bar
    %     c.Label.FontSize = 20;                                              % font size of label of color bar
    %     fig_name = sprintf('Re_stress%d',k);
    %     saveas(gcf,fig_name ,'png'); 
    %     % saveas(gcf,'fig' ,'png'); 
    % end  
    % skewness(歪度)を算出
    % if k < 15
    % if k == 6 || k == 7 || k == 8|| k == 14
    %     figure(n)
    %     grid on
    %     grid minor
    %     box on
    %     xlabel('\sly^{+}','FontName','Times','FontAngle','Italic','FontSize',20);
    %     ylabel('\slS','FontName','Times','FontAngle','Italic','FontSize',20);
    %     xlim([0 60])
    %     set(gca,'FontName','Times','FontSize',15)
    %     hold on
    %     ycross = x(1:n,1,2) * ut / vis;
    %     sigma(:,:,1)      = sqrt(uu(:,:,1));
    %     sigma(:,:,2)      = sqrt(uu(:,:,2));
    %     skewness_u = (uu_Skewness(:,:,1)./(sigma(:,:,1).^3));
    %     skewness_v = (uu_Skewness(:,:,2)./(sigma(:,:,2).^3));
    %     % skewness_uの第二引数はどのx座標を観察したいかによって適宜変える
    %     semilogx(ycross,skewness_u(:,1),'--');
    %     % ylim([0 60])
    %     hold off
    % end 
    % flatnessを算出
    % if k < 15
    % if k == 6 || k == 7 || k == 8|| k == 14
    %     figure(n)
    %     grid on
    %     grid minor
    %     box on
    %     xlabel('\sly^{+}','FontName','Times','FontAngle','Italic','FontSize',20);
    %     ylabel('\slF','FontName','Times','FontAngle','Italic','FontSize',20);
    %     xlim([0 60])
    %     set(gca,'FontName','Times','FontSize',15)
    %     hold on
    %     ycross = x(1:n,1,2) * ut / vis;
    %     sigma(:,:,1)      = sqrt(uu(:,:,1));
    %     sigma(:,:,2)      = sqrt(uu(:,:,2));
    %     flatness_u = (uu_Flatness(:,:,1)./(sigma(:,:,1).^4));
    %     flatness_v = (uu_Flatness(:,:,2)./(sigma(:,:,2).^4));

    %     % skewness_uの第二引数はどのx座標を観察したいかによって適宜変える
    %     semilogx(ycross,flatness_u(:,1),'o');
    %     % ylim([0 60])
    %     hold off
    % end     
end

    uus(:,1)      = sqrt(mean(uu(:,:,1),2));
    uus(:,2)      = sqrt(mean(uu(:,:,2),2));
    uus(:,3)      = mean(uu(:,:,3),2);
    euv(1:n-1,1)  = uus(1:n-1,3)./diff(ums(:,1)).*diff(x(:,1,2));
% 
% 1D-figure
figure;
subplot(2,2,1);
% p = plot(x(1:n,1,2)/h,uus(:,1)/ut,'ko');
p = plot(x(1:n,1,2)*ut/vis,uus(:,1)/ut,'ko');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.1 1 1];
box on;
%axis equal tight;                                                 
% xticks([0:0.2:1]);
% yticks([0:0.5:2.5]);
% xlim([0 1]);
xlim([0 60]); 
% ylim([0 2.5]);
ylim([0 5.0]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
% xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
ylabel('$${\it u^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');
% 
subplot(2,2,2);
solution_plot_Xaxis = x(1:n,1,2)*ut/vis
solution_u_RMS_plot_Yaxis = uus(:,1)/ut
solution_v_RMS_plot_Yaxis = uus(:,2)/ut
% RMSのグラフ作成に必要な変数を保存
save('session_solutionRMS','solution_plot_Xaxis','solution_u_RMS_plot_Yaxis', 'solution_v_RMS_plot_Yaxis')
p = plot(x(1:n,1,2)*ut/vis,uus(:,2)/ut,'k^');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [1 0.2 1];
box on;
%axis equal tight;                                                 
% xticks([0:0.2:1]);
% yticks([0:0.5:1.5]);
% xlim([0 1]);
% ylim([0 1.5]);
xlim([0 60]); 
ylim([0 5.0]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
% xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
ylabel('$${\it v^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');

% レイノルズせん断応力
subplot(2,2,3);
solution_Normalized_Re_stress = uus(:,3)/(ut*ut)
save('solution_Normalized_Re_stress','solution_Normalized_Re_stress')
% p = plot(x(1:n,1,2)/h,uus(:,3)/(ut*ut),'k>');
p = plot(x(1:n,1,2)*ut/vis,uus(:,3)/(ut*ut),'k>');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.5 0.5 1];
box on;
%axis equal tight;                                                 
% xticks([0:0.2:1]);
% yticks([-1:0.5:1.0]);
% xlim([0 1]);
% ylim([0 1]);
xlim([0 60]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
% xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
ylabel('$$-\overline{{\it u^\prime}^+{\it v^\prime}^+}$$','FontSize',20,'Interpreter','latex');
% 
subplot(2,2,4);
p = plot(x(1:n-1,1,2)/h,euv(:,1)/vis,'ko');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.1 0.1 1];
box on;
%axis equal tight;                                                 
xticks([0:0.2:1]);
yticks([0:20:80]);
xlim([0 1]);
ylim([0 80]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
ylabel('$$\epsilon_{u^\prime v^\prime}^+$$','FontSize',20,'Interpreter','latex');
saveas(gcf,'Fluc-y','png');
% 
% %
    toc;

