A = u_time_space_sum;
B = fixed_coord_y(:,1)*10^-3;
figure;plot(A,B);
scatter(A,B);
title('velocity distribution of mean mainstream direction')
xlabel('u[m/s]')
ylabel('y[m]')
