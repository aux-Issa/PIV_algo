close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
tic;

% variables
m        = 63;                      % number of y direction
n        = 53;                      % number of y direction
num_data = 100;                     % number of velocity-files

ro = 997.05;                        %   density
h = 0.02;                           %   half chanel width[m]
U_b = 0.335;                        %   Bulk velocity[m/s]

path_of_coordinate=('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/preTest_Re15000.6v7z21dj/analyzed_data/coordinate/coordinate.6v7z8lhi.000000.dat');
coord = importdata(path_of_coordinate);
xy = coord.data; % separating text and data in struct
xy = xy-xy(1,:); % set x1 to zero
y(:,:,1) = flip(col2im(xy(:,1),[m n],[m n],'distinct'),2); % x data [mm] with flip
y(:,:,2) = col2im(xy(:,2),[m n],[m n],'distinct'); % y data [mm]

% x座標
y = y/1000; 
% y座標
y = y/1000;  

% velocity_data: 速度データを一つのファイルに格納
velocity_data = cell(1,num_data);
for j = 0:num_data-1
   if j<10
       file_name = sprintf('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/preTest_Re15000.6v7z21dj/analyzed_data/velocity/velocity_Re15000.6v7z8lhi.00000%d.dat',j);
   elseif j<100
       file_name = sprintf('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/preTest_Re15000.6v7z21dj/analyzed_data/velocity/velocity_Re15000.6v7z8lhi.0000%d.dat',j);
   elseif j<1000
       file_name = sprintf('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/preTest_Re15000.6v7z21dj/analyzed_data/velocity/velocity_Re15000.6v7z8lhi.000%d.dat',j);
   end      
   raw_velocity_data{1,j+1} = readmatrix(file_name);
end

for i = 1:num_data
    velocity_data{1,i}(:,1) = rmmissing(raw_velocity_data{1,i}(:,1));
    velocity_data{1,i}(:,2) = rmmissing(raw_velocity_data{1,i}(:,2));
end
velocity_data
disp('velocity_data') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 速度uの時間平均
time_u_sum = 0;
for i=1:num_data
    time_u_sum = time_u_sum+velocity_data{1,i}(:,1);
end    
raw_time_u_mean = time_u_sum / num_data;

% 速度uの時間平均を座標系に整形
time_u_mean = zeros(n,m);
initialIndex = 1;
for i=1:n
 time_u_mean(i , :) = raw_time_u_mean(initialIndex : (initialIndex + m - 1) , 1);
 initialIndex = initialIndex + m;
end    

u_time_space_sum = 0;
for i=1:m
    u_time_space_sum = u_time_space_sum + time_u_mean(:,i);
end 
u_time_space_sum = u_time_space_sum/m;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = u_time_space_sum/U_b;
A = y(1,1:n,2)/h;
figure;plot(A,B);
scatter(A,B);
title('velocity distribution of mean mainstream direction')

xlabel('y/h')
ylabel('u/U_b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 速度vの時空間平均
time_v_sum = 0;
for i=1:num_data
%     time_v_sum = time_v_sum+fixed_velocity_data{1,i}(:,1);
time_v_sum = time_v_sum+velocity_data{1,i}(:,2);
end    
time_v_mean=time_v_sum/num_data;
% OK
exact_time_v_mean = zeros(n,m);
initialIndex = 1;
for i=1:n
 exact_time_v_mean(i , :) = time_v_mean(initialIndex : (initialIndex + m - 1) , 1);
 initialIndex = initialIndex + m;
 initialIndex;
end    
% fixed_time_v_mean = reshape(time_v_mean,[m,n]);
% fixed_coord_y = reshape(coord.data(:,2),[m,n]);
% 
v_time_space_sum = 0;
for i=1:n
    v_time_space_sum = v_time_space_sum + exact_time_v_mean(i, :);
end 
v_time_space_average = v_time_space_sum/n;

B_v = v_time_space_average/U_b;
A_v = y(1:m,1,1)/h;
figure;plot(A_v,B_v);
scatter(A_v,B_v);
title('velocity distribution of mean mainstream direction')

xlabel('y/h')
ylabel('u/U_b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%速度uのデータのみを取り出し
u_data = zeros(m*n,num_data);
for i= 1:num_data
    u_data(:,i) = rmmissing(velocity_data{1,i}(:,1));
end   


u_dash = zeros(m*n,num_data);
for i=1:num_data
    for v= 1:m*n
        u_dash(v,i) = u_dash(v,i) + velocity_data{1,i}(v,1) - raw_time_u_mean(v,1);
    end
end 
v_dash = zeros(m*n,num_data);

for i=1:num_data
    for v= 1:m*n
        v_dash(v,i) = u_dash(v,i) + velocity_data{1,i}(v,2)-time_v_mean(v,1);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%レイノルズせん断応力の式の右辺、u_dashとv_dashの積uv_dashを算出
uv_dash = u_dash.*v_dash;
uv_sum = 0;
for i=1:num_data
    uv_sum = uv_sum + uv_dash(:,i);
end
%uv_dashの時間平均を算出
uv_dash_mean = uv_sum/num_data; 
% fixed_uv_dash_mean = reshape(uv_dash_mean, [m,n]);

uv_dash_sum = zeros(n,m);
initialIndex = 1;
for i=1:n
    uv_dash_sum(i , :) = uv_dash_mean(initialIndex : initialIndex + m  -1, 1);
    initialIndex = initialIndex + m;
end
%x軸のデータ数で除してuv_dashの時空間平均を算出
space_uv_dash_sum = 0;
for i=1:m
    space_uv_dash_sum = space_uv_dash_sum + uv_dash_sum(:,i);
end
space_uv_dash_mean = space_uv_dash_sum/m;

%レイノルズせん断応力 τ=-ρ*space_uv_dash_meanを算出する
re_stress = -(ro*10^-9)*space_uv_dash_mean;
% re_stress = -(ro*10^-9)*space_uv_dash_mean;

A = y(1,1:n,2)/h;
B = re_stress;
% B = -space_uv_dash_mean/(U_b^2);
 figure;plot(A,B(:));

 