clear all
close all
clc

tic
%Aluminum
% rho = 2710; %kg/m^3
% c = 1256;    %J/kg-C
% k = 180;    %W/m-K

% %Copper
% rho = 8954;
% k = 386;
% c = 380;

%Carbon Steel
% rho = 7753;
% k = 36;
% c = 486;

%Gold
% rho = 18900;
% k = 318;
% c = 130;

%Silver
rho = 10510;
k = 418;
c = 230;
        
del = 0.0025;%meters
dt = 0.99*(c*rho*del^2*del^2)/(4*k*del^2);
tfinal = 10; %seconds
tvec = 0:dt:tfinal;
n_times = length(tvec);

%Create Main Shape
[x,y] = meshgrid(0:del:0.1,0:del:0.2);
[n_y,n_x] = size(x);
temp_one = ones(n_y,n_x);
temp_zero = zeros(n_y,n_x);
temp_shape1 = y > -sqrt(0.05^2 - (x-0.05).^2)+0.05;
temp_shape2 = y < sqrt(0.05^2 - (x-0.05).^2)+0.15;
temp_shape = temp_one .* (temp_shape1 & temp_shape2);

%Create the constant temperature condition
temp_side = temp_zero;
temp_side(:,[1,end]) = temp_shape(:,[1,end]);

figure
mesh(x,y,temp_side)


%Create the round boundary conditions
temp_curve = temp_zero;
temp_outside = temp_zero;
for col = 2:n_x-1
    slice_ind = find(temp_shape(:,col));
    temp_curve([max(slice_ind),min(slice_ind)],col) = 1;
end

for row = 1:n_y
    slice_ind = find(temp_shape(row,:));
    temp_curve(row,[max(slice_ind),min(slice_ind)]) = 1;
end
temp_curve = temp_curve - temp_side;

temp_normal = temp_shape-temp_side-temp_curve;
figure
mesh(x,y,temp_normal)

figure
mesh(x,y,temp_curve)

%Plot the geometry
figure
mesh(x,y,temp_shape)


%%
%Initialize u
u = zeros(n_y,n_x,n_times);
u(:,:,1) = temp_shape*60;
C = (k*dt/(c*rho*del^2));

[row_side,col_side] = find(temp_side);

[row_normal,col_normal] = find(temp_normal);

[row_curve,col_curve] = find(temp_curve);
I_movie = 1;
flux = 100;

for I_time = 1:n_times-1
    for I_x = 1:n_x
        for I_y = 1:n_y
            %Constant Temp Condition
            if 1 == sum(I_y == row_side & I_x == col_side)
                u(I_y,I_x,I_time+1) = 80;
            
            %Curved Edge Condition
            elseif 1 == sum(I_y == row_curve & I_x == col_curve)
                if I_y > n_y/2
                    n = [x(I_y,I_x)-0.05,y(I_y,I_x)-0.15];
                else
                    n = [x(I_y,I_x)-0.05,y(I_y,I_x)-0.05];
                end
                n = n./norm(n);
                
                if (n(1) >= 0) && (n(2) >= 0) %Top right
                    ux = u(I_y,I_x-1,I_time)-2*flux*n(1)*del; %u_x+1
                    uy = u(I_y-1,I_x,I_time)-2*flux*n(2)*del; %u_y+1
                    u(I_y,I_x,I_time+1) = u(I_y,I_x,I_time)+C*(ux-2*u(I_y,I_x,I_time)+u(I_y,I_x-1,I_time))+C*(uy-2*u(I_y,I_x,I_time)+u(I_y-1,I_x,I_time));
                elseif (n(1) < 0) && (n(2) >= 0) %Top left
                    ux = u(I_y,I_x+1,I_time)+2*flux*n(1)*del; %u_x-1
                    uy = u(I_y-1,I_x,I_time)-2*flux*n(2)*del; %u_y+1
                    u(I_y,I_x,I_time+1) = u(I_y,I_x,I_time)+C*(u(I_y,I_x+1,I_time)-2*u(I_y,I_x,I_time)+ux)+C*(uy-2*u(I_y,I_x,I_time)+u(I_y-1,I_x,I_time));
                elseif (n(1) < 0) && (n(2) < 0) %Bottom left
                    ux = u(I_y,I_x+1,I_time)+2*flux*n(1)*del; %u_x-1
                    uy = u(I_y+1,I_x,I_time)+2*flux*n(2)*del; %u_y-1
                    u(I_y,I_x,I_time+1) = u(I_y,I_x,I_time)+C*(u(I_y,I_x+1,I_time)-2*u(I_y,I_x,I_time)+ux)+C*(u(I_y+1,I_x,I_time)-2*u(I_y,I_x,I_time)+uy);
                else %Bottom Right
                    ux = u(I_y,I_x-1,I_time)-2*flux*n(1)*del; %u_x+1
                    uy = u(I_y+1,I_x,I_time)+2*flux*n(2)*del; %u_y-1
                    u(I_y,I_x,I_time+1) = u(I_y,I_x,I_time)+C*(ux-2*u(I_y,I_x,I_time)+u(I_y,I_x-1,I_time))+C*(u(I_y+1,I_x,I_time)-2*u(I_y,I_x,I_time)+uy);
                end
              
            %Inside Condition
            elseif 1 == sum(I_y == row_normal & I_x == col_normal)
                u(I_y,I_x,I_time+1) = u(I_y,I_x,I_time)+C*(u(I_y,I_x+1,I_time)-2*u(I_y,I_x,I_time)+u(I_y,I_x-1,I_time))+C*(u(I_y+1,I_x,I_time)-2*u(I_y,I_x,I_time)+u(I_y-1,I_x,I_time));
            else
            end
            
            
        end
    end
    if 2 == mod(I_time,10)
                
                fig = figure(5);
                surf(x,y,u(:,:,I_time))
                M(I_movie) = getframe(fig)
                I_movie = I_movie+1;
    else
    end
end

toc

%movie2avi(M, 'aluminum.avi')
%%
figure(6)
surf(x,y,u(:,:,round(0.5/dt)))

figure(7)
surf(x,y,u(:,:,round(2/dt)))

figure(8)
surf(x,y,u(:,:,round(5/dt)))

figure(9)
surf(x,y,u(:,:,round(10/dt)))

