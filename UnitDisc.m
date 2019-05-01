L = 1;
T = 30;
H = 0.1;
rho = 1;
alpha = H/rho;
dt = 0.9 * dx * sqrt(rho/H);
dx = 0.01;
dy = 0.01;
PosX = unique([0:dx:L,L]);
Time = unique([0:dt:T,T]);
PosY = unique([0:dy:L,L]);
Sol = zeros(length(PosX), length(PosY), length(Time));
for k = 2: length(PosY)
    Sol(:,k,2) = (PosX.^2 + PosY.^2)/10;
%     Sol(:,k,2) = (exp(-10*(PosX-(L/2)).^2) - exp(-10*(L/2).^2))*dt;
end
r = alpha * dt^2;
index = 2:length(PosX)-1;
ind = 2:length(PosY)-1;

for l = index
    a = PosX(l);
    b = PosY(l);
    if ((a^2 + b^2) == 1)
        Sol(l, l, :) = 0.1;
    end
    if ((a^2 + b^2) > 1)
        Sol(l, l, :) = 0;
    end
end

for j = 3: length(Time) 
    Sol(index, ind, j) = 2*Sol(index, ind, j-1) - Sol(index, ind, j-2) + r/dx^2*(Sol(index+1, ind, j-1)) - 2*Sol(index, ind, j-1) + Sol(index-1, ind, j-1) + r/dy^2*(Sol(index, ind+1, j-1)) - 2*Sol(index, ind, j-1) + Sol(index, ind-1, j-1);
    for l = index
        a = PosX(l);
        b = PosY(l);
        if ((a^2 + b^2) == 1)
            Sol(l, l, :) = 0.1;
        end
        if ((a^2 + b^2) > 1)
            Sol(l, l, :) = 0;
        end
    end
end

[x, timeframes] = size(Time);
figure(1);
for i = 1:3
    h = surf(PosX, PosY, Sol(:,:,i)');
    set(h, 'LineStyle', 'none');
    hold on;
end
hold off
title('Wave over a two-dimension plate');
xlabel('x(position)'); ylabel('y(position)'); zlabel('u(displacement)');