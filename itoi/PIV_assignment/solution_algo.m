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

    uf(:,:,1)  = u(:,:,1)-ums(:,1);
    uf(:,:,2)  = u(:,:,2)-ums(:,2);
    uf(:,:,3)  = -uf(:,:,1).*uf(:,:,2);

    uu(:,:,1)  = uu(:,:,1)+uf(:,:,1).*uf(:,:,1)/double(numfiles);
    uu(:,:,2)  = uu(:,:,2)+uf(:,:,2).*uf(:,:,2)/double(numfiles);
    uu(:,:,3)  = uu(:,:,3)+uf(:,:,3)/double(numfiles);
    U_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,1),xq,yq,'spline');         % interpolation
    V_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,2),xq,yq,'spline');         % interpolation

    if k < 15
        [cc_solution_dash,hc_solution_dash]=contourf(xq*ut/vis,yq*ut/vis,U_dash/ub,8);                   % making isoline
        hc_solution_dash.TextStep = 0.4;                                                  % interval isoline
        hc_solution_dash.ShowText = 'off';                                                % isoline text off
        colormap('jet');                                                    % color type of 'jet'
        box on;                                                             % making flame of figure
        xlim([0 150]);                                                       % range of x
        ylim([0 60]);                                                     % range of y
        set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
        ax = gca;
        ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
        xlabel('${\it x^+}$','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
        ylabel('${\it y^+}$','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
        c = colorbar;                                                       % making color bar
        c.Limits = [-0.4 0.4];                                                 % range of colorbar
        c.FontSize = 18;                                                    % font size of scale of color bar
        c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
        c.Label.Interpreter = 'latex';                                      % font type of label of color bar
        c.Label.String = '$${\it u^\prime}/{\it U_b}$$';                % label of color bar
        c.Label.FontSize = 20;                                              % font size of label of color bar
        fig_name = sprintf('U_dash_solution%d',k);
        saveas(gcf,fig_name ,'png'); 
        
    end 
    if k < 15
        [cc_solution_dash,hc_solution_dash]=contourf(xq*ut/vis,yq*ut/vis,V_dash/ub,8);                   % making isoline
        hc_solution_dash.TextStep = 0.4;                                                  % interval isoline
        hc_solution_dash.ShowText = 'off';                                                % isoline text off
        colormap('jet');                                                    % color type of 'jet'
        box on;                                                             % making flame of figure
        xlim([0 150]);                                                       % range of x
        ylim([0 60]);                                                     % range of y
        set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
        ax = gca;
        ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
        xlabel('${\it x^+}$','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
        ylabel('${\it y^+}$','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
        c = colorbar;                                                       % making color bar
        c.Limits = [-0.4 0.4];                                                 % range of colorbar
        c.FontSize = 18;                                                    % font size of scale of color bar
        c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
        c.Label.Interpreter = 'latex';                                      % font type of label of color bar
        c.Label.String = '$${\it u^\prime}/{\it U_b}$$';                % label of color bar
        c.Label.FontSize = 20;                                              % font size of label of color bar
        fig_name = sprintf('U_dash_solution%d',k);
        saveas(gcf,fig_name ,'png'); 
        
    end  
end

    uus(:,1)      = sqrt(mean(uu(:,:,1),2));
    uus(:,2)      = sqrt(mean(uu(:,:,2),2));
    uus(:,3)      = mean(uu(:,:,3),2);
    euv(1:n-1,1)  = uus(1:n-1,3)./diff(ums(:,1)).*diff(x(:,1,2));
% 
% 1D-figure
figure;
subplot(2,2,1);
p = plot(x(1:n,1,2)/h,uus(:,1)/ut,'ko');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.1 1 1];
box on;
%axis equal tight;                                                 
xticks([0:0.2:1]);
yticks([0:0.5:2.5]);
xlim([0 1]);
ylim([0 2.5]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
ylabel('$${\it u^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');
% 
subplot(2,2,2);
p = plot(x(1:n,1,2)/h,uus(:,2)/ut,'k^');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [1 0.2 1];
box on;
%axis equal tight;                                                 
xticks([0:0.2:1]);
yticks([0:0.5:1.5]);
xlim([0 1]);
ylim([0 1.5]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
ylabel('$${\it v^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');

% レイノルズせん断応力
subplot(2,2,3);
p = plot(x(1:n,1,2)/h,uus(:,3)/(ut*ut),'k>');
p.LineWidth = 1.2;
p.MarkerSize = 7;
p.MarkerEdgeColor = 'black';
p.MarkerFaceColor = [0.5 0.5 1];
box on;
%axis equal tight;                                                 
xticks([0:0.2:1]);
yticks([-1:0.5:1.0]);
xlim([0 1]);
ylim([0 1]);
set( gca, 'FontName','Times','FontSize',18); 
ax = gca;
ax.TickLength = [0.02 0.1];
xlabel('{\it y}/{\it h}','FontSize',20,'Interpreter','latex');
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

