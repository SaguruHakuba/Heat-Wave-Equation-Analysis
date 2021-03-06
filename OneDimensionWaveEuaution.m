L = 2;
T = 30;
H = 0.1;
rho = 1;
alpha = H/rho;
dt = 0.9 * dx * sqrt(rho/H);
dx = 0.01;
PosX = unique([0:dx:L,L]);
Time = unique([0:dt:T,T]);
Sol = zeros(length(PosX), length(Time));
Sol(:,2) = (exp(-10*(PosX-(L/2)).^2) - exp(-10*(L/2).^2))*dt;
r = alpha * dt^2 / dx^2;
index = 2:length(PosX)-1;
for j = 3: length(Time)
    Sol(index, j) = 2*Sol(index, j-1) - Sol(index, j-2) + r*(Sol(index+1, j-1)) - 2*Sol(index, j-1) + Sol(index-1, j-1);
end

[x, timeframes] = size(Time);
figure(1);
for i = 1:10
    plot(PosX, Sol(:,i)');
    hold on;
end
hold off
title('One dimension');
xlabel('x(position)'); ylabel('u(displacement)');