% % 2D-figure
figure;
newpoints = m*3;
% % クエリグリッドを作成
% [xq,yq] = meshgrid(...                                              % making fine mesh
%                linspace(x(1,1,1),x(1,m,1),newpoints ),...
%                linspace(x(1,1,2),x(n,1,2),newpoints )...
%              );
%  クエリ点で内挿
U_dash= interp2(x(:,:,1),x(:,:,2),uf(:,:,1),xq,yq,'spline');         % interpolation
[cc_water_dash,hc_water_dash]=contourf(xq*ut/vis,yq*ut/vis,U_dash/ut,20);                   % making isoline
hc_water_dash.TextStep = 0.4;                                                  % interval isoline
hc_water_dash.ShowText = 'off';                                                % isoline text off
colormap('jet');                                                    % color type of 'jet'
box on;                                                             % making flame of figure
xlim([0 300]);                                                       % range of x
ylim([0 100]);                                                     % range of y
set( gca, 'FontName','Times','FontSize',18);                        % font and its size of axes 
ax = gca;
ax.TickLength = [0.02 0.1];                                         % scale size to inside from flame
xlabel('{\it x}^+','FontSize',20,'Interpreter','latex');           % xlabel, its size, and type
ylabel('{\it y}^+','FontSize',20,'Interpreter','latex');           % ylabel, its size, and type
c = colorbar;                                                       % making color bar
c.Limits = [0 3.0];                                                 % range of colorbar
c.FontSize = 18;                                                    % font size of scale of color bar
c.TickLabelInterpreter = 'latex';                                   % font type of scale of color bar
c.Label.Interpreter = 'latex';                                      % font type of label of color bar
c.Label.String = '$${\it \u_dash}/{\it \overline{u^+}}$$';                % label of color bar
c.Label.FontSize = 20;                                              % font size of label of color bar
saveas(gcf,'U_dash','png');   