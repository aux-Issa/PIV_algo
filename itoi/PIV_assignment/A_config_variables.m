close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
tic;
m        = 63;    
n        = 53;     % number of y direction
num_data = 100;

ro = 997.05; %density
h = 0.02; %half chanel width[m]
U_b = 0.335; %Bulk velocity[m/s]

loadname=('/Volumes/HDCZ-UT/itoi_PIV/water/water_test.6uvaasgh/preTest_Re15000.6v7z21dj/analyzed_data/coordinate/coordinate.6v7z8lhi.000000.dat');
coord = importdata(loadname);
xy = coord.data; % separating text and data in struct
xy = xy-xy(1,:); % set x1 to zero
x(:,:,1) = flip(col2im(xy(:,1),[m n],[m n],'distinct'),2); % x data [mm] with flip
x(:,:,2) = col2im(xy(:,2),[m n],[m n],'distinct'); % y data [mm]
x = x/1000; 
velocity_data=cell(1,num_data);
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
%速度uを正に変換
for i = 1:num_data
    velocity_data{1,i}(:,1) = rmmissing(raw_velocity_data{1,i}(:,1));
    velocity_data{1,i}(:,2) = rmmissing(raw_velocity_data{1,i}(:,2));
end
