water_plot_Xaxis = load('session_waterRMS', 'water_plot_Xaxis').water_plot_Xaxis;
water_u_RMS_plot_Yaxis = load('session_waterRMS', 'water_u_RMS_plot_Yaxis').water_u_RMS_plot_Yaxis;
water_v_RMS_plot_Yaxis = load('session_waterRMS', 'water_v_RMS_plot_Yaxis').water_v_RMS_plot_Yaxis;
solution_plot_Xaxis = load('session_solutionRMS', 'solution_plot_Xaxis').solution_plot_Xaxis;
solution_u_RMS_plot_Yaxis = load('session_solutionRMS', 'solution_u_RMS_plot_Yaxis').solution_u_RMS_plot_Yaxis;
solution_v_RMS_plot_Yaxis = load('session_solutionRMS', 'solution_v_RMS_plot_Yaxis').solution_v_RMS_plot_Yaxis;
water_Normalized_Re_stress = load('water_Normalized_Re_stress', 'water_Normalized_Re_stress').water_Normalized_Re_stress;
solution_Normalized_Re_stress = load('solution_Normalized_Re_stress', 'solution_Normalized_Re_stress').solution_Normalized_Re_stress;

% p_u_water = plot(water_plot_Xaxis,water_u_RMS_plot_Yaxis,'b^');
% p_u_water.LineWidth = 1.2;
% p_u_water.MarkerSize = 7;
% % p_u_water.MarkerEdgeColor = 'blue';
% p_u_water.MarkerFaceColor = 'blue';
% box on;
% xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
% ylabel('$${\it u^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');
% xlim([0 60]); 
% hold on
% p_u_solution = plot(solution_plot_Xaxis,solution_u_RMS_plot_Yaxis,'r^');
% p_u_solution.MarkerEdgeColor = 'red';
% p_u_solution.MarkerFaceColor = 'red';
% legend({'water','80ppm'},'Location','northeast')
% hold off 

% p_v_water = plot(water_plot_Xaxis,water_v_RMS_plot_Yaxis,'b^');
% p_v_water.LineWidth = 1.2;
% p_v_water.MarkerSize = 7;
% p_v_water.MarkerFaceColor = 'blue';
% box on;
% xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
% ylabel('$${\it v^\prime}_{{\rm rms}}^+$$','FontSize',20,'Interpreter','latex');
% xlim([0 60]); 
% hold on
% p_v_solution = plot(solution_plot_Xaxis,solution_v_RMS_plot_Yaxis,'r^');
% p_v_solution.MarkerEdgeColor = 'red';
% p_v_solution.MarkerFaceColor = 'red';
% legend({'water','80ppm'},'Location','northeast')
% hold off 

% 無次元化されたレイノルズ剪断応力の粘弾性底層内の比較
p_v_water = plot(water_plot_Xaxis,water_Normalized_Re_stress,'b^');
p_v_water.LineWidth = 1.2;
p_v_water.MarkerSize = 7;
p_v_water.MarkerFaceColor = 'blue';
box on;
xlabel('$${\it y^+}$$','FontSize',20,'Interpreter','latex');
ylabel('$$-\overline{{\it u^\prime}^+{\it v^\prime}^+}$$','FontSize',20,'Interpreter','latex');
xlim([0 60]); 
hold on
p_v_solution = plot(solution_plot_Xaxis,solution_Normalized_Re_stress,'r^');
p_v_solution.MarkerEdgeColor = 'red';
legend({'water','80ppm'},'Location','northeast')
hold off 