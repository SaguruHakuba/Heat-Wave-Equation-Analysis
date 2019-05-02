% Brandon M. Keltz CM3229
% Project 4
% MA342 - Computational Modeling
% April 29, 2019
% I pledge this is my work.

dx = 0.01;
dy = 0.01;
r = 0.1;
dt = r * min(dx, dy);
T = 5;
p = inf;
alpha = 2;

x = -1:dx:1;
y = 1:-dy:-1;
[X, Y] = meshgrid(x, y);

if p == inf
	
	pu = max(abs(X), abs(Y)) / 10;
	pu(max(abs(X), abs(Y)) >= 1) = 0.1;
	
else
	
	pu = ((abs(X).^p + abs(Y).^p).^(1 / p)) / 10;
	pu(abs(X).^p + abs(Y).^p >= 1) = 0.1;
	
end

cu = pu;
nu = cu;

for t = 0:dt:T
	
% 	figure(1);
% 	if p == inf
% 		strtitle = ['Solution at time, \textit{t} = ', num2str(t, '%.3f'), ', where \textit{p} = $\infty$'];
% 	else
% 		strtitle = ['Solution at time, \textit{t} = ', num2str(t, '%.3f'), ', where \textit{p} = ', num2str(p, '%d')];
% 	end
% 	s = surf(X, Y, cu);
% 	set(s, 'LineStyle', 'none');
% 	title(strtitle, 'Interpreter', 'latex');
% 	xlabel('\textit{x}', 'Interpreter', 'latex');
% 	ylabel('\textit{y}', 'Interpreter', 'latex');
% 	zlabel('\textit{u(x, y)}', 'Interpreter', 'latex');
% 	axis([-1, 1, -1, 1, 0, 0.2]);
% 	pause(0.000000000000000001);
	
	nu(2:end - 1, 2:end - 1) = 2 * (1 - alpha * dt^2 / dx^2 - alpha * dt^2 / dy^2) * cu(2:end - 1, 2:end - 1) ...
		+ alpha * dt^2 * ((cu(2:end - 1, 1:end - 2) + cu(2:end - 1, 3:end)) / dx^2 ...
		+ (cu(1:end - 2, 2:end - 1) + cu(3:end, 2:end - 1)) / dy^2) - pu(2:end - 1, 2:end - 1);
	
	if p == inf
		
		nu(max(abs(X), abs(Y)) >= 1) = 0.1;
		
	else
		
		nu(X.^p + Y.^p >= 1) = 0.1;
		
	end
	
	pu = cu;
	cu = nu;
	
end

s = surf(X, Y, cu);
set(s, 'LineStyle', 'none');
xlabel('\textit{x}', 'Interpreter', 'latex');
ylabel('\textit{y}', 'Interpreter', 'latex');
zlabel('\textit{u(x, y)}', 'Interpreter', 'latex');
axis([-1, 1, -1, 1, 0, 0.2]);