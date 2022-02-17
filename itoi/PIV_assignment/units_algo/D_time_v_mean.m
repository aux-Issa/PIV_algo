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
 initialIndex
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
A_v = x(1:m,1,1)/h;
figure;plot(A_v,B_v);
scatter(A_v,B_v);
title('velocity distribution of mean mainstream direction')

xlabel('y/h')
ylabel('u/U_b')


