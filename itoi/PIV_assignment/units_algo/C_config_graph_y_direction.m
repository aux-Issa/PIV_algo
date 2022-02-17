B = u_time_space_sum/U_b;
A = x(1,1:n,2)/h;
figure;plot(A,B);
scatter(A,B);
title('velocity distribution of mean mainstream direction')

xlabel('y/h')
ylabel('u/U_b')

