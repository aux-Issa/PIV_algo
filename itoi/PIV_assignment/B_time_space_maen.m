time_u_sum = 0;
for i=1:num_data
%     time_u_sum = time_u_sum+fixed_velocity_data{1,i}(:,1);
time_u_sum = time_u_sum+velocity_data{1,i}(:,1);
end    
time_u_mean=time_u_sum/num_data;
% OK
exact_time_u_mean = zeros(n,m);
initialIndex = 1;
for i=1:n
 exact_time_u_mean(i , :) = time_u_mean(initialIndex : (initialIndex + m - 1) , 1);
 initialIndex = initialIndex + m;
 initialIndex
end    
% fixed_time_u_mean = reshape(time_u_mean,[m,n]);
% fixed_coord_y = reshape(coord.data(:,2),[m,n]);
% 
time_space_sum = 0;
 for i=1:m
     time_space_sum = time_space_sum + exact_time_u_mean(:,i);
 end 
 time_space_mean = time_space_sum/m;

