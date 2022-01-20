

u_dash_squre = u_dash.* u_dash;
u_dash_squre_sum = 0;

for i= 1:num_data
    u_dash_squre_sum = u_dash_squre_sum + u_dash_squre(:,i);
end
u_dash_squre_mean = u_dash_squre_sum/500;
fixed_u_dash_squre_mean = reshape(u_dash_squre_mean,[128,253]);
fixed_u_dash_squre_mean_sqrt = sqrt(fixed_u_dash_squre_mean);

sum=0;
for i=1:n
    sum = sum + fixed_u_dash_squre_mean_sqrt(:,i);
end
rms_u_dash = sum/n;
A = rms_u_dash;
B = fixed_coord_y*10^-3;
figure;plot(A,B);
scatter(A,B);
title('rms-distribution of velocity fluctation in the x-direction')
xlabel('rms of u-velocity fluctation[m/s]')
ylabel('y[m]')

%--------------------------------------------------------------------------
v_dash_squre = v_dash.* v_dash;
v_dash_squre_sum = 0;

for i= 1:num_data
    v_dash_squre_sum = v_dash_squre_sum + v_dash_squre(:,i);
end
v_dash_squre_mean = v_dash_squre_sum/500;
fixed_v_dash_squre_mean = reshape(v_dash_squre_mean,[128,253]);
fixed_v_dash_squre_mean_sqrt = sqrt(fixed_v_dash_squre_mean);
sum=0;
for i=1:n
    sum = sum + fixed_v_dash_squre_mean_sqrt(:,i);
end
rms_v_dash = sum/n;
fixed_coord_x = reshape(coord.data(:,1),[128,253]);

A = rms_v_dash;
B = fixed_coord_y(:,1)*10^-3;
figure;plot(A,B);
scatter(A,B);
title('rms-distribution of velocity fluctation in the y-direction')
xlabel('rms of v-velocity fluctation[m/s]')
ylabel('y[m]')