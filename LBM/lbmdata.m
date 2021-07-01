function [U, V, compare] = lbmdata(nx,Re,x,y,u,ux,uy)
  
% compare with Ghia's paper (reference)
% by kimhaemulgae
  
if nx == 129
    u_velocity = flipud(ux([1 8 9 10 14 23 37 59 65 80 95 110 123 124 125 126 129],65))./u;
    v_velocity =flipud(uy(65,[1 9 10 11 13 21 30 31 65 104 111 117 122 123 124 125 129])')./u;
elseif nx == 257
    u_velocity = flipud(ux([1 15 17 19 27 45 73 117 129 159 189 219 245 247 249 251 257],129))./u;
    v_velocity =flipud(uy(129,[1 17 19 21 25 41 59 61 129 207 221 233 243 245 247 249 257])')./u;
elseif nx == 513
    u_velocity = flipud(ux([1,29,33,37,53,89,145,233,257,317,377,437,489,493,497,501,513],257))./u;
    v_velocity =flipud(uy(257,[1,33,37,41,49,81,117,121,257,413,441,465,485,489,493,497,513])')./u;
elseif nx == 1025
    u_velocity = flipud(ux([1,57,65,73,105,177,289,465,513,633,753,873,977,985,993,1001,1025],513))./u;
    v_velocity =flipud(uy(513,[1,65,73,81,97,161,233,241,513,825,881,929,969,977,985,993,1025])')./u;

end

if Re == 100
    u_thesis = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.21090 -0.15662 -0.10150 -0.06434 -0.04775 -0.04192 -0.03717 0]';
    v_thesis = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.10890 0.10091 0.09233 0]';
elseif Re == 1000
    u_thesis = [1 0.65928 0.57492 0.51117 0.46604 0.33304 0.18719 0.05702 -0.06080 -0.10648 -0.27805 -0.38289 -0.29730 -0.22220 -0.20196 -0.18109 0]';
    v_thesis = [0 -0.21388 -0.27669 -0.33714 -0.39188 -0.51550 -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 0.30353 0.29012 0.27485 0]';
% elseif Re == 10
%     u_thesis = 
end

error_u = abs((u_thesis - u_velocity)./u_thesis)*100;
error_v = abs((v_thesis - v_velocity)./v_thesis)*100;

if Re==100
    error_u_y_04531 = error_u(10)
    error_v_x_08047 = error_v(8)
    error_v_x_02344 = error_v(10)
    uy_min = u_velocity(10);
    vx_min = v_velocity(8);
    vx_max = v_velocity(10);
    compare = [uy_min; vx_min; vx_max] 
elseif Re==1000
    error_u_y_01719 = error_u(12);
    error_v_x_09063 = error_v(6);
    error_v_x_01563 = error_v(12);
    uy_min = u_velocity(12)
    vx_min = v_velocity(6)
    vx_max = v_velocity(12)
    compare = [uy_min; vx_min; vx_max] 
elseif Re==10
    
end 

y_axis = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0]';
x_axis = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0]';

% figure(2),
% % plot(y_axis,u_thesis,y_axis,u_thesis,'o'),title('U'), grid on, hold on, plot(y_axis, u_velocity,'r'), legend('Paper', 'Paper', 'Calculated'); 
% plot(flipud(y(:,1)),flipud(ux(:,(nx+1)/2)./u),'r'), xlabel('y','fontsize',15), ylabel('ux','fontsize',15), hold on,plot(y_axis,u_thesis,'o','MarkerFaceColor','b'), legend('LBM', 'Ghia','Location','Best');
% figure(3),
% % plot(x_axis,v_thesis,x_axis,v_thesis,'o'),title('V'), grid on, hold on, plot(x_axis, v_velocity,'r'), legend('Paper', 'Paper', 'Calculated'); 
% plot(flipud(x(1,:)),flipud(uy((nx+1)/2,:)./u),'r'), xlabel('x','fontsize',15), ylabel('uy','fontsize',15), hold on,plot(x_axis,v_thesis,'o','MarkerFaceColor','b'), legend('LBM', 'Ghia','Location','Best');

U = flipud(ux(:,(nx+1)/2)./u);
V = flipud(uy((nx+1)/2,:)./u)';

end




