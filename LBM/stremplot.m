function streamplot(NN,ux,uy)

% streamline plotting  
  
    [x,y] = meshgrid(0:1/(NN-1):1,0:1/(NN-1):1);

    [startx,starty] = meshgrid(0:0.1:1,0:0.1:1);

%     starty=.49:.005:.51;
%     startx=.49*ones(length(starty),1);
% 
%     starty=.475:.002:.525;
%     startx=.475*ones(length(starty),1);

    streamline(x,y,ux,uy,startx,starty);

    xlim([0 1])
    ylim([0 1])

%     xlim([0.37 0.49])
%     ylim([.47 .53])
end
