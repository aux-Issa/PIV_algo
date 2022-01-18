time_u_sum = 0;
for i=1:num_data
    time_u_sum = time_u_sum+fixed_velocity_data{1,i}(:,1);
end    
time_u_mean=time_u_sum/500;
fixed_time_u_mean = reshape(time_u_mean,[128,253]);
fixed_coord_y = reshape(coord.data(:,2),[128,253]);

time_space_sum = 0;
for i=1:253
    time_space_sum = time_space_sum + fixed_time_u_mean(:,i);
end 
time_space_mean = time_space_sum/253;

