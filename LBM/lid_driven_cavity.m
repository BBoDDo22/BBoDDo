%% Lattice Boltzmann Method (LBM) _ 2D lid_driven_cavity by kimhaemulgae

clear all, close all, clc;

NN = 257;
xmin = 0; xmax = 1; nx = NN; x = linspace(xmin, xmax, nx);
dx = x(2)-x(1);
ymin = 0; ymax = 1; ny = NN; y = linspace(ymin, ymax, ny);
dy = y(2)-y(1); 

% unit lattice condition
c=1; dt=dx/c;
Re = 100; 
u = 0.1; 
nu = u*256/Re;
tau = 3*nu + 0.5;

% c=1; dt = dx/c; dt=1; Re = 100;
% tau = 1/1.25; nu = (tau-0.5)/3;
% u = nu*Re/(nx-1);

[x, y] = meshgrid(x, y);

nodenums = nx*ny; % lattice node numbers
left  = [1:nx:nodenums]';
right = [nx:nx:nodenums]';
bottom = [1:1:nx]';
top = [nx*(ny-1)+1:1:nodenums]';

% lattice speed c
% c=1; dt=dx/c; cs=c;
% c = dx/dt;
% cs = c/sqrt(3);

% particle velocity vectors e
ec = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1];
e = c*ec;
eT = e';

% inflow velocity and density
inflow_u = ones(size(top))*[u 0]; 
dens = 2.7;

% initial condition
U = zeros(nodenums, 2);
U(top,:) = inflow_u; 
U(bottom,:) = 0;
U(right,:) = 0;
U(left,:) = 0; 


% equilibrium distribution function of particles
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
opp = [1 4 5 2 3 8 9 6 7]; 
eU = (U*eT);
U2 = sum(U.^2, 2);
ro = dens*ones(nodenums, 1);
f_initial = ones(nodenums, 1) * w;
feq_initial = (ro*w).*(1 + 3*eU/c^2 + (9/2)*(eU.^2)/c^4 - (3/2)*U2*ones(1,9)/(c^2)); 
fstar = feq_initial;

% relaxation parameter 
% tau = 1; nu = ((2*tau-1)*dx^2)/(dt*6); Re = (inflow_u(1,1)*(ymax-ymin)*1)/nu
% nu = 10; tau = (3*nu)/(c^2*dt) + 1/2; Re = inflow_u(1,1)*(ymax-ymin)/nu

% shift arrays
Right = [2:nx]; Left = [1:nx-1]; 
Top = [2:ny];  Bottom = [1:ny-1]; 

u_convergence = 1;
v_convergence = 1;
ttt = 1;
tic
while u_convergence >= 1e-9 || v_convergence >= 1e-9
% for z=1:300000
    
    % step3
    ro  = sum(fstar,2);
    U = (fstar*e)./(ro*ones(1,2));
    
    U(top,:) = inflow_u;
    U(left,:) = 0; U(bottom,:) = 0;
    U(right,:) = 0;
    ro(top) = (1./(1+U(top,2))).*(fstar(top,1)+fstar(top,2)+fstar(top,4)+2*(fstar(top,3)+fstar(top,6)+fstar(top,7)));

    % Zou-He
    fstar(top,5) = fstar(top,3) - 2/3*ro(top).*U(top,2);
    fstar(top,8) = fstar(top,6) + 1/2*(fstar(top,2)-fstar(top,4)) - 1/2*ro(top).*U(top,1);
    fstar(top,9) = fstar(top,7) - 1/2*(fstar(top,2)-fstar(top,4)) + 1/2*ro(top).*U(top,1);

    % step4
    eU = (U*eT);
    U2 = sum(U.^2, 2);
    feq = (ro*w).*(1 + 3*eU/c^2 + (9/2)*(eU.^2)/c^4 - (3/2)*U2*ones(1,9)/(c^2));

    % step5 : collision
    f = fstar - (1/tau)*(fstar - feq);
        
%     % save boundary for bounceback
%     fright = fstar(right, [2 6 9]);
%     ftop = fstar(top,[3 6 7]);    
%     fbottom = fstar(bottom, [5 8 9]);

    % mid grid BC
    f(right,[4 8 7]) = fstar(right, [2 6 9]);
    f(bottom,[3 6 7]) = fstar(bottom, [5 8 9]);
    f(left,[2 6 9]) = fstar(left, [4 8 7]);
    
%     f(right,:) = fstar(right, opp);
%     f(top,:) = fstar(top, opp);
%     f(bottom,:) = fstar(bottom, opp);


    
    % Stream
    [fstar] = Stream(fstar, f, nodenums, nx, ny, Top, Bottom, Right, Left);    
    
%     % mid grid BC
%     fstar(right,[4 8 7]) = f(right, [2 6 9]);
%     fstar(bottom,[3 6 7]) = f(bottom, [5 8 9]);
%     fstar(left,[2 6 9]) = f(left, [4 8 7]);
% 
   
%     % Zou-He
%     fstar(top,5) = fstar(top,3) - 2/3*ro(top).*U(top,2);
%     fstar(top,8) = fstar(top,6) + 1/2*(fstar(top,2)-fstar(top,4)) - 1/2*ro(top).*U(top,1);
%     fstar(top,9) = fstar(top,7) - 1/2*(fstar(top,2)-fstar(top,4)) + 1/2*ro(top).*U(top,1);
 
%     fstar(top,3) = fstar(top,5);
%     fstar(top,6) = fstar(top,8) + 1/2*(fstar(top,4)-fstar(top,2)) + 1/2*ro(top).*U(top,1);
%     fstar(top,7) = fstar(top,9) + 1/2*(fstar(top,2)-fstar(top,4)) - 1/2*ro(top).*U(top,1);

    
    ux = reshape(U(:,1), nx, ny)';
    uy = reshape(U(:,2), nx, ny)';
    V = sqrt(ux.^2 + uy.^2);
    
%     if tic == 10
%         preV = V;
%     end
%     
%     if tic>10
%     error = sum(sum(V-preV))/sum(sum(V));
%     preV = V;
%     end

    if ttt == 10
        pre_ux = ux;
        pre_uy = uy;
    end
    
%     if ttt>10
%         convergence = sum(sum((ux-pre_ux).^2 + (uy-pre_uy).^2))/sum(sum(pre_ux.^2+pre_uy.^2));
%         convergence = sqrt(convergence);
%         pre_ux = ux;
%         pre_uy = uy;
%     end

    if ttt>10
        u_convergence = abs(sum(sum((ux-pre_ux).^2))^0.5/sum(sum(pre_ux)));
        v_convergence = abs(sum(sum((uy-pre_uy).^2))^0.5/sum(sum(pre_uy)));
        pre_ux = ux;
        pre_uy = uy;
    end
        
    
    if(mod(ttt, 100) == 0) 
% %         subplot(4,2,[3:8]), 
%         figure(1), imagesc(V); axis equal;set(gca,'Ydir','Normal'); title('Velocity Re100 257','fontsize',15); colorbar;
% %         figure(2), quiver(ux,uy), axis equal
% %         subplot(421); imagesc(ux); axis equal; title('x-velocity','fontsize',15);
% %         subplot(422); imagesc(uy); axis equal; title('y-velocity','fontsize',15);        
        figure(1),subplot(121), imagesc(V); axis equal;set(gca,'Ydir','Normal'); title('Velocity Re100 257','fontsize',15); colorbar; colormap(jet);
        Ro = reshape(ro,NN,NN)';
        subplot(122),imagesc(Ro); axis equal;set(gca,'Ydir','Normal'); title('Ro Re100 257','fontsize',15); colorbar; colormap(jet);

        drawnow;

        ttt
        u_convergence, v_convergence
    end
    ttt=ttt+1;
     
    
end
toc
CalculatedTime=toc;
CalculatedHour=CalculatedTime/3600
[U,V,compare] = lbmdata(nx,Re,x,y,u,ux,uy);
% save('cavity_Re100_257')

streamplot(NN,ux,uy)


% lx = linspace(xmin, xmax, nx/5);
% ly = linspace(xmin, xmax, nx/5);
% % figure(2), quiver(ux,uy,1,'b-')
% figure, streamline(x,y,ux,uy,lx,ly)
% figure, streamline(x(1:129,1:129),y(1:129,1:129),ux(1:129,1:129),uy(1:129,1:129),lx,ly)
% figure, streamline(x(1:129,129:end),y(1:129,129:end),ux(1:129,129:end),uy(1:129,129:end),lx,ly)

