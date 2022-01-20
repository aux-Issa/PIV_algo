%速度uのデータのみを取り出し
u_data = zeros(m*n,num_data);
for i= 1:num_data
    u_data(:,i) = rmmissing(velocity_data{1,i}(:,1));
end   


u_dash = zeros(m*n,num_data);
for i=1:num_data
    for v= 1:m*n
        u_dash(v,i) = u_dash(v,i) + velocity_data{1,i}(v,1)-time_u_mean(v,1);
    end
end 
v_dash = zeros(m*n,num_data);

for i=1:num_data
    for v= 1:m*n
        v_dash(v,i) = u_dash(v,i) + velocity_data{1,i}(v,2)-time_v_mean(v,1);
    end
end 
