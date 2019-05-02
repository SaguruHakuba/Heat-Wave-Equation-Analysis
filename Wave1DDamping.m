% Brandon M. Keltz CM3229
% Project 4
% MA342 - Computational Modeling
% April 29, 2019
% I pledge this is my work.

dx = 0.01;
L = 2;
H = 0:0.01:1;
rho = 1;
r = 0.1;
dt = 0.1 * dx;
T = 1:100;
alpha = 0;
beta = 0;
epsilon = 1e-3;
x = (0:dx:L)';

kappa = zeros(length(T), length(H));
% hold on;

j = 1;

for h = H
	
	for i = 1:length(T)
		
		akappa = 0;
		bkappa = 10;
		
		while abs(bkappa - akappa) > epsilon
			
			mkappa = (bkappa + akappa) / 2;
			
			pu = zeros(length(x), 1);
			cu = pu + dt * (exp(-10 * (x - L / 2).^2) - exp(-10 * (L / 2).^2));
			
			cu(1) = alpha;
			cu(end) = beta;
			nu = cu;
			
			t = 0;
			
			while all(nu >= 0) && t <= T(i)
				
				%         figure(1);
				%         strtitle = ['Solution at time, \textit{t} = ', num2str(i, '%.3f'), ', where $\kappa = $ ', num2str(mkappa, '%.3f')];
				%         plot(x, cu);
				%         title(strtitle, 'Interpreter', 'latex');
				%         xlabel('\textit{x}', 'Interpreter', 'latex');
				%         ylabel('\textit{y}', 'Interpreter', 'latex');
				%         axis([0, L, 0, 1]);
				%         pause(0.001);
				nu(2:end - 1) = ((2 - 2 * H(j) * r^2 / rho + mkappa * dt) * cu(2:end - 1) + H(j) * r^2 * (cu(1:end - 2) + cu(3:end)) / rho - pu(2:end - 1)) / (1 + mkappa * dt);
				nu(1) = alpha;
				nu(end) = beta;
				
				pu = cu;
				cu = nu;
				t = t + dt;
				
			end
			
			%         plot(i, bkappa, 'r*');
			%         plot(i, akappa, 'b*');
			
			if any(nu < 0)
				
				akappa = mkappa;
				
			else
				
				bkappa = mkappa;
				
			end
			
		end
		
		%     hold off;
		kappa(i, j) = mkappa;
		
	end
	
	j = j + 1;
	
end