close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
tic;
m        = 128;    
n        = 253;     % number of y direction
num_data = 500;

ro = 997.05; %density


loadname=('/Users/issa/Desktop/HeatTransfer/assignment/入門課題2/COORD.dat');
coord = importdata(loadname);
xy = coord.data; % separating text and data in struct
xy = xy-xy(1,:); % set x1 to zero
x(:,:,1) = flip(col2im(xy(:,1),[m n],[m n],'distinct'),2); % x data [mm] with flip
x(:,:,2) = col2im(xy(:,2),[m n],[m n],'distinct'); % y data [mm]
x = x/1000; 
data=cell(1,num_data);
for j = 0:num_data-1
   if j<10
       file_name = sprintf('/Users/issa/Desktop/HeatTransfer/assignment/入門課題2/DATA/water.somehow1.00000%d.dat',j);
   elseif j<100
       file_name = sprintf('/Users/issa/Desktop/HeatTransfer/assignment/入門課題2/DATA/water.somehow1.0000%d.dat',j);
   elseif j<1000
       file_name = sprintf('/Users/issa/Desktop/HeatTransfer/assignment/入門課題2/DATA/water.somehow1.000%d.dat',j);
   end      
   data{1,j+1} = readmatrix(file_name);
end
%速度uを正に変換
for i = 1:num_data
   fixed_velocity_data{1,i}(:,1) = rmmissing(data{1,i}(:,1)*(-1));
   fixed_velocity_data{1,i}(:,2) = rmmissing(data{1,i}(:,2));
end
