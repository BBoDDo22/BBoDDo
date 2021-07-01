%%  uniform flow by kimhaemulgae

clear all, close all, clc;

% required function : ForceDistributionFunc, IBM, Stream
NN = 801;
xmin = 0; xmax = 1; nx = NN*xmax - (xmax-1); x = gpuArray(linspace(xmin, xmax, nx));
dx = 1/(NN-1);
ymin = 0; ymax = 1; ny = NN; y = gpuArray(linspace(ymin, ymax, ny));
dy = 1/(NN-1); 
 
phy_u = 1e-3; % m/s
phy_l = 1e-3; % m
phy_D = 1e-3/40; % m
phy_density = 1000; % kg/m^3
phy_dx = phy_l/(NN-1); % physic's length/lattice Numbers
Lag_r = phy_D/2/phy_l; r=Lag_r; d=2*r; % 
Lag_dens = phy_density*phy_l^3;
phy_nu = 1e-3;

Lattice_Unit_dx = 1; Lattice_Unit_dt = 1; Lattice_Unit_u = 0.1; Lattice_Unit_dens = 1;
Cx = phy_dx/Lattice_Unit_dx;
Cu = phy_u/Lattice_Unit_u; Ct = Cx/Cu; 
Cdens = phy_density/Lattice_Unit_dens;

% Dimension corrections
phy_dt = Lattice_Unit_dt*Ct;
Lattice_Unit_D = phy_D/Cx;
Lattice_Unit_r = Lattice_Unit_D/2;

Re_cylinder = 40,   Title = ['Cylinder Re' num2str(Re_cylinder) ' ' num2str(NN)];
Lattice_Unit_nu = Lattice_Unit_u*2*(Lattice_Unit_r)/Re_cylinder;
Lattice_Unit_tau = 3*Lattice_Unit_nu + 0.5
Re=Lattice_Unit_u*2*nx/Lattice_Unit_nu;

c = Lattice_Unit_dx/Lattice_Unit_dt;
in_u = Lattice_Unit_u; 
nu = Lattice_Unit_nu;
tau = Lattice_Unit_tau;
dens = Lattice_Unit_dens;
dt = Lattice_Unit_dt;

[x, y] = meshgrid(x, y);

nodenums = nx*ny; % lattice node numbers
left  = gpuArray([1:nx:nodenums]');
right = gpuArray([nx:nx:nodenums]');
bottom = gpuArray([1+1:1:nx-1]');
top = gpuArray([nx*(ny-1)+1+1:1:nodenums-1]');

% particle velocity vectors e
e = gpuArray([0 0; c 0; 0 c; -c 0; 0 -c; c c; -c c; -c -c; c -c]);
eT = e';

% inflow velocity and density
inflow_u = gpuArray(ones(size(left))*[in_u 0]); 

% initial condition
U = gpuArray(zeros(nodenums, 2));
F = gpuArray(zeros(nodenums,9)); Eux = gpuArray(zeros(nx,ny)'); Euy = Eux;
desired_velocity = 0;
U (:,1) = Lattice_Unit_u;

% equilibrium distribution function of particles
w = gpuArray([4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]);
W = ones(nodenums,1)*w; % nodenums x 9
opp = gpuArray([1 4 5 2 3 8 9 6 7]); 
eU = (U*eT)/c;
U2 = sum(U.^2, 2);
ro = gpuArray(dens*ones(nodenums, 1));
f_initial = ones(nodenums, 1) * w;
feq = (ro*w).*(1 + 3*eU/c^2 + (9/2)*(eU.^2)/c^4 - (3/2)*U2*ones(1,9)/(c^2)); 
fstar = feq;

% shift arrays
Right = gpuArray([2:nx]); Left = gpuArray([1:nx-1]); 
Top = gpuArray([2:ny]);  Bottom = gpuArray([1:ny-1]); 

% IBM
% initial Lagrangian points:
dd=2*pi*Lattice_Unit_r/(1/1.5);
Ldx = gpuArray((0:1/dd:1)');
circle = [r*cos(2*pi*Ldx(1:end))+0.45 r*sin(2*pi*Ldx(1:end))+ymax/2]; % point positions
Lb_nodenums = length(circle); % circle boundary node 개수
L_center = sum(circle)/Lb_nodenums; % center
Larea = (2*pi*(Lattice_Unit_r)/(Lb_nodenums));
 
% initial fluid flow, force:
u0  = gpuArray(zeros(Lb_nodenums,2)); 

% declare fluid variables
Lux = gpuArray(zeros(Lb_nodenums,1));
Luy = Lux;
u  = u0;

% Lagrangian points:
Lx = gpuArray(circle(:,1)); Ly = gpuArray(circle(:,2));
Lfx = gpuArray(zeros(Lb_nodenums,1)); Lfy = gpuArray(zeros(Lb_nodenums,1));
Efx = gpuArray(zeros(nx,ny)); Efy = gpuArray(zeros(nx,ny)); 
fib = gpuArray(zeros(nodenums, 2));
Eux = gpuArray(zeros(nx,ny)); Euy = gpuArray(zeros(nx,ny));
L = circle; circ = zeros(nx,ny); circc = zeros(nx,ny); Circ = zeros(Lb_nodenums,1);
desired_velocity = zeros(Lb_nodenums,2);

%% first forcing step

error = 1;
boundary_error = 1;
nt = 500000;
ttt=1;
videotitle = [Title '_test1.avi'];
% v = VideoWriter(videotitle); open(v);
% parpool(4);
% spmd
tic
while error >= 1e-6
% for z=1:100
%     z
    
    % step4 : collision  
    F = ForceDistributionFunc(nodenums,e,U,fib,tau,w);
    f = fstar - (1/tau)*(fstar - feq) + F*Lattice_Unit_dt;    
     
%     % bounce back
%     f(top,[5 8 9]) = f(top, [3 6 7]);
%     f(bottom,[3 6 7]) = f(bottom, [5 8 9]);

    % step1 Stream
    fstar = Stream(fstar, f, nodenums, nx, ny, Top, Bottom, Right, Left);
      
    % Zou-He B.C.           
    U(left,:) = inflow_u; % u=1, v=0
%     U(left,1) = poise'; U(left,2)=0; % u=1, v=0
    ro(left) = (1./(1-U(left,1))).*(fstar(left,1)+fstar(left,3)+fstar(left,5)+2*(fstar(left,4)+fstar(left,7)+fstar(left,8)));    
    fstar(left,2) = fstar(left,4) + 2/3*ro(left).*U(left,1);
    fstar(left,6) = fstar(left,8) - 1/2*(fstar(left,3)-fstar(left,5)) + 1/6*ro(left).*U(left,1) + 1/2*ro(left).*U(left,2);
    fstar(left,9) = fstar(left,7) + 1/2*(fstar(left,3)-fstar(left,5)) + 1/6*ro(left).*U(left,1) - 1/2*ro(left).*U(left,2);
    
    U(right,:) = 2*U(right-1,:)-U(right-2,:); % dv/dx=0
%     U(right,:) = (4*U(right-1,:)-U(right-2,:))/3; % dv/dx=0
    fstar(right,2) = 2*fstar(right-1,2)-fstar(right-2,2); % second
    fstar(right,6) = 2*fstar(right-1,6)-fstar(right-2,6);
    fstar(right,9) = 2*fstar(right-1,9)-fstar(right-2,9);
%     fstar(right,2) = (4*fstar(right-1,2)-fstar(right-2,2))/3; % third
%     fstar(right,6) = (4*fstar(right-1,6)-fstar(right-2,6))/3;
%     fstar(right,9) = (4*fstar(right-1,9)-fstar(right-2,9))/3;

    
%     U(top,1) = 0; 
%     U(top,1) = (2*U(top-nx,1)-U(top-2*nx,1)); % du/dy=0
    U(top,1) = (4*U(top-nx,1)-U(top-2*nx,1))/3; % du/dy=0
    U(top,2) = 0; % v=0   
%     U(top,:) = (4*U(top-nx,:)-U(top-2*nx,:))/3;
    ro(top) = (1./(1+U(top,2))).*(fstar(top,1)+fstar(top,2)+fstar(top,4)+2*(fstar(top,3)+fstar(top,6)+fstar(top,7)));
    fstar(top,5) = fstar(top,3) - 2/3*ro(top).*U(top,2);
    fstar(top,8) = fstar(top,6) + 1/2*(fstar(top,2)-fstar(top,4)) - 1/2*ro(top).*U(top,1);
    fstar(top,9) = fstar(top,7) - 1/2*(fstar(top,2)-fstar(top,4)) + 1/2*ro(top).*U(top,1);
% 
    U(bottom,1) = 0;
%     U(bottom,1) = (2*U(bottom+nx,1)-U(bottom+2*nx,1)); % du/dy=0
    U(bottom,1) = (4*U(bottom+nx,1)-U(bottom+2*nx,1))/3; % du/dy=0
    U(bottom,2) = 0; % v=0
%     U(bottom,:) = (4*U(bottom+nx,:)-U(bottom+2*nx,:))/3;
    ro(bottom) = (1./(1-U(bottom,2))).*(fstar(bottom,1)+fstar(bottom,2)+fstar(bottom,4)+2*(fstar(bottom,5)+fstar(bottom,8)+fstar(bottom,9)));
    fstar(bottom,3) = fstar(bottom,5) + 2/3*ro(bottom).*U(bottom,2);
    fstar(bottom,6) = fstar(bottom,8) - 1/2*(fstar(bottom,2)-fstar(bottom,4)) + 1/2*ro(bottom).*U(bottom,1);
    fstar(bottom,7) = fstar(bottom,9) + 1/2*(fstar(bottom,2)-fstar(bottom,4)) - 1/2*ro(bottom).*U(bottom,1);

%     % Periodic
%     fstar(bottom,[3 6 7]) = fstar(top, [3 6 7]);
%     fstar(top,[5 8 9]) = fstar(bottom, [5 8 9]);

%     % modified bounce back
%     fstar(top,[5 8 9]) = f(top, [3 6 7]);
%     fstar(bottom,[3 6 7]) = f(bottom, [5 8 9]);
%     fstar(top,:) = f(top, opp);
%     fstar(bottom,:) = f(bottom, opp);

    ro  = sum(fstar,2);  
    U = (fstar*e)./(ro*ones(1,2));

    % Immersed Boundary Method Lu, Lf, Ef    
    Eux = reshape(U(:,1), nx, ny);
    Euy = reshape(U(:,2), nx, ny);
       
    Lux = gpuArray(zeros(Lb_nodenums,1));
    Luy = Lux;
    Efx = gpuArray(zeros(nx,ny)); Efy = gpuArray(zeros(nx,ny));
    
    [Lux Luy Lfx Lfy Efx Efy fib desired_velocity old_Lux old_Luy R] = IBM_GPU2(Lx, Ly, desired_velocity, Eux, Euy, Lux, Lfx, Efx, Luy, Lfy, Efy, fib, ro, dx, dy, dt, Larea, nodenums,nx, ny, Lb_nodenums,Lattice_Unit_dx);
%     [Lux Luy Lfx Lfy Efx Efy fib desired_velocity old_Lux old_Luy R] = MDF_IBM_GPU(Lx, Ly, desired_velocity, Eux, Euy, Lux, Lfx, Efx, Luy, Lfy, Efy, fib, ro, dx, dy, dt, Larea, nodenums,nx, ny, Lb_nodenums,Lattice_Unit_dx,U);
   
    U = U + (fib*Lattice_Unit_dt)./(2*ro*ones(1,2));
    
    eU = (U*eT);
    U2 = sum(U.^2, 2);
    feq = (ro*w).*(1 + 3*eU/c^2 + (9/2)*(eU.^2)/c^4 - (3/2)*U2*ones(1,9)/(c^2));

    old_Lux = zeros(Lb_nodenums,1);
    old_Luy = zeros(Lb_nodenums,1);
           
    Eux = reshape(U(:,1), nx, ny);
    Euy = reshape(U(:,2), nx, ny);
    V = sqrt(Eux.^2 + Euy.^2)';
    
    cir = find(Efx~=0);

    CdE(ttt) = abs(sum(fib(:,1)))/(mean(ro)*Lattice_Unit_u^2*(Lattice_Unit_r));
    CdL(ttt) = abs(sum(fib(:,2)))/(mean(ro)*Lattice_Unit_u^2*(Lattice_Unit_r));

    if(mod(ttt, 100) == 0)
        ttt
        CD = [CdE(ttt), CdL(ttt)]
        
        figure(1), 
        imagesc(V); axis equal;set(gca,'Ydir','Normal'); colorbar;colormap(jet); title([Title '  CdE ' num2str(CD(1)) ' CdL ' num2str(CD(2))],'fontsize',10); colorbar;colormap(jet);        
        
%         frame = getframe(gcf);
%         writeVideo(v,frame);
% 
%         drawnow;
        
        if ttt == 10000
            pre_ux = Eux(10:end-10,10:end-10);
            pre_uy = Euy(10:end-10,10:end-10);
        end
        if ttt>10000
            error = sum(sum((Eux(10:end-10,10:end-10)-pre_ux).^2 + (Euy(10:end-10,10:end-10)-pre_uy).^2))/sum(sum(pre_ux.^2+pre_uy.^2));
            error = sqrt(error)
            pre_ux = Eux(10:end-10,10:end-10);
            pre_uy = Euy(10:end-10,10:end-10);
        end
                
        if error >= 1e-3 && error <= 1.1e-3
            e3error=error;
            e3Cd = CD;
            e3ttt = ttt;
        end

        if error >= 1e-4 && error <= 1.1e-4
            e4error=error;
            e4Cd = CD;
            e4ttt= ttt;
        end
        
        if error >= 1e-5 && error <= 1.1e-5
            e5error=error;
            e5Cd = CD;
            e5ttt = ttt;
        end

        if ttt >= 10000
            boundary_error = sqrt(sum((Lux.^2 + Luy.^2))/Lb_nodenums)
        end        

    end
              
    
    ttt=ttt+1;
    
end

e6error=error;
e6Cd = CD;
e6ttt = ttt;

toc
