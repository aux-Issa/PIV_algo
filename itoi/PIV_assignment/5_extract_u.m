%速度uのデータのみを取り出し
u_data = zeros(32384,500);
for i= 1:num_data
    u_data(:,i) = rmmissing(data{1,i}(:,1));
end   


u_dash = zeros(3284,500);
for i=1:num_data
    for v= 1:m*n
        u_dash(:,i) = u_dash + data{1,i}(v,1)-time_u_mean(v,1)
    end
end 