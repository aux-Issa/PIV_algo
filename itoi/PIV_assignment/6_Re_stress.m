%レイノルズせん断応力の式の右辺、u_dashとv_dashの積uv_dashを算出
uv_dash = u_dash.*v_dash;
uv_sum = 0;
for i=1:500
    uv_sum = uv_sum + uv_dash(:,i);
end
%uv_dashの時間平均を算出
uv_dash_mean = uv_sum/500; 
fixed_uv_dash_mean = reshape(uv_dash_mean, [128,253]);
space_uv_dash_sum = 0;
for i=1:n
    space_uv_dash_sum = space_uv_dash_sum + fixed_uv_dash_mean(:,1);
end
%x軸のデータ数で除してuv_dashの時空間平均を算出
space_uv_dash_mean = space_uv_dash_sum/n;


%レイノルズせん断応力 τ=-ρ*space_uv_dash_meanを算出する
re_stress = -(ro*10^-9)*space_uv_dash_mean
%plot
A = re_stress;
B = fixed_coord_y;
figure;plot(A,B(:));