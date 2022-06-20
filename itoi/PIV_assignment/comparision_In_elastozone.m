water_plot_Xaxis = load('session_waterRMS', 'water_plot_Xaxis').water_plot_Xaxis;
water_u_RMS_plot_Yaxis = load('session_waterRMS', 'water_u_RMS_plot_Yaxis').water_u_RMS_plot_Yaxis;
water_v_RMS_plot_Yaxis = load('session_waterRMS', 'water_v_RMS_plot_Yaxis').water_v_RMS_plot_Yaxis;
solution_plot_Xaxis = load('session_solutionRMS', 'solution_plot_Xaxis').solution_plot_Xaxis;
solution_u_RMS_plot_Yaxis = load('session_solutionRMS', 'solution_u_RMS_plot_Yaxis').solution_u_RMS_plot_Yaxis;
solution_v_RMS_plot_Yaxis = load('session_solutionRMS', 'solution_v_RMS_plot_Yaxis').solution_v_RMS_plot_Yaxis;
water_Normalized_Re_stress = load('water_Normalized_Re_stress', 'water_Normalized_Re_stress').water_Normalized_Re_stress;
solution_Normalized_Re_stress = load('solution_Normalized_Re_stress', 'solution_Normalized_Re_stress').solution_Normalized_Re_stress;
DNSfile_name = sprintf('/Users/issa/Desktop/DNS.csv');
DNSdata = csvread(DNSfile_name);
y_cross_DNS = DNSdata(:,2);
u_mean_DNS = DNSdata(:,3);
u_RMS_DNS = DNSdata(:,4).^(1/2);
v_RMS_DNS = DNSdata(:,5).^(1/2);
f = 1;
% DNS_lowRe_file_name = sprintf('/Users/issa/Desktop/DNS_Re.csv');
% DNS_lowRe_data = csvread(DNS_lowRe_file_name);
% y_cross_lowRe_DNS = DNS_lowRe_data(:,2);
% u_mean_lowRe_DNS = DNS_lowRe_data(:,3);
% u_RMS_lowRe_DNS = DNS_lowRe_data(:,4).^(1/2);
% v_RMS_lowRe_DNS = DNS_lowRe_data(:,5).^(1/2);

%% Uのrms
figure(f)
p_u_water = plot(water_plot_Xaxis,water_u_RMS_plot_Yaxis,'b^');
p_u_water.LineWidth = 1.2;
p_u_water.MarkerSize = 7;
p_u_water.MarkerFaceColor = 'blue';
box on;
xlabel('$${\it y^+}$$','FontSize',30,'Interpreter','latex');
ylabel('$${\it u^\prime}_{{\rm rms}}^+$$','FontSize',30,'Interpreter','latex');
xlim([0 60]); 
ylim([0 5]);
hold on
p_u_solution = plot(solution_plot_Xaxis,solution_u_RMS_plot_Yaxis,'ro');
p_u_solution.MarkerEdgeColor = 'red';
p_u_solution.MarkerFaceColor = 'red';
p_u_DNS = plot(y_cross_DNS,u_RMS_DNS,'k-');
% p_u_lowRe_DNS = plot(y_cross_lowRe_DNS,u_RMS_lowRe_DNS,'k^');

legend({'water','80ppm','water-DNS'},'Location','northeast')
hold off 
f = f + 1;

%% Vのrms
figure(f)
p_v_water = plot(water_plot_Xaxis,water_v_RMS_plot_Yaxis,'b^');
p_v_water.LineWidth = 1.2;
p_v_water.MarkerSize = 7;
p_v_water.MarkerFaceColor = 'blue';
box on;
xlabel('$${\it y^+}$$','FontSize',30,'Interpreter','latex');
ylabel('$${\it v^\prime}_{{\rm rms}}^+$$','FontSize',30,'Interpreter','latex');
xlim([0 60]); 
ylim([0 5]);
hold on
p_v_solution = plot(solution_plot_Xaxis,solution_v_RMS_plot_Yaxis,'ro');
p_v_solution.MarkerEdgeColor = 'red';
p_v_solution.MarkerFaceColor = 'red';
p_v_DNS = plot(y_cross_DNS,v_RMS_DNS,'k-');
legend({'water','80ppm'},'Location','northeast')
hold off 
f = f + 1;

% 無次元化されたレイノルズ剪断応力の粘弾性底層内の比較
figure(f)
p_v_water = plot(water_plot_Xaxis,water_Normalized_Re_stress,'b^');
p_v_water.LineWidth = 1.2;
p_v_water.MarkerSize = 7;
p_v_water.MarkerEdgeColor = 'blue';
p_v_water.MarkerFaceColor = 'blue';
box on;
xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
ylabel('$$-\overline{{\it u^\prime}^+{\it v^\prime}^+}$$','FontSize',20,'Interpreter','latex');
xlim([0 60]); 
ylim([-0.6 1]); 
yticks([-1:0.2:1.0]);
hold on
p_v_solution = plot(solution_plot_Xaxis,solution_Normalized_Re_stress,'ro');
p_v_solution.MarkerEdgeColor = 'red';
p_v_solution.MarkerFaceColor = 'red';
legend({'water','80ppm'},'Location','northeast')
hold off 
toc;
% todo: 関数定義 and (他の言語に置き換え AWS S3とかにCSVアップロードできてもいいかも)