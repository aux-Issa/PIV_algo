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
space_uv_dash_mean = space_uv_dash_sum / m;

%レイノルズせん断応力 τ=-ρ*space_uv_dash_meanを算出する
re_stress = -(ro*10^-9)　*　space_uv_dash_mean;
% re_stress = -(ro*10^-9)*space_uv_dash_mean;

A = x(1,1:n,2)　/　h;
B = re_stress;
% B = -space_uv_dash_mean/(U_b^2);
 figure;plot(A,B(:));

 