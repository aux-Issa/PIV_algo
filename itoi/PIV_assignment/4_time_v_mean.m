%速度vの時間平均
time_v_sum = 0;
for i=1:num_data
    time_v_sum = time_v_sum + fixed_velocity_data{1,i}(:,2);
end   
time_v_mean = time_v_sum/500;
fixed_time_v_mean = reshape(time_v_mean,[128,253]);

   