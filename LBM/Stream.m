function [fstar] = Stream(fstar, f, nodenums, nx, ny, Top, Bottom, Right, Left)

% streaming step function by kimhaemulgae

for i=1:9
        if i==1
            fstar1=reshape(fstar(:,i),nx,ny)'; f1=reshape(f(:,i),nx,ny)';
            fstar1=f1;
            fstar1=reshape(fstar1',nodenums,1);
        elseif i==2
            fstar2=reshape(fstar(:,i),nx,ny)'; f2=reshape(f(:,i),nx,ny)';
            fstar2(:,Right) = f2(:,Left);
            fstar2=reshape(fstar2',nodenums,1);
        elseif i==3
            fstar3=reshape(fstar(:,i),nx,ny)'; f3=reshape(f(:,i),nx,ny)';
            fstar3(Top,:) = f3(Bottom,:);
            fstar3=reshape(fstar3',nodenums,1);
        elseif i==4
            fstar4=reshape(fstar(:,i),nx,ny)'; f4=reshape(f(:,i),nx,ny)';
            fstar4(:,Left) = f4(:,Right);
            fstar4=reshape(fstar4',nodenums,1);
        elseif i==5
            fstar5=reshape(fstar(:,i),nx,ny)'; f5=reshape(f(:,i),nx,ny)';
            fstar5(Bottom,:) = f5(Top,:);
            fstar5=reshape(fstar5',nodenums,1);
        elseif i==6
            fstar6=reshape(fstar(:,i),nx,ny)'; f6=reshape(f(:,i),nx,ny)';
            fstar6(Top,Right) = f6(Bottom,Left);
            fstar6=reshape(fstar6',nodenums,1);
        elseif i==7
            fstar7=reshape(fstar(:,i),nx,ny)'; f7=reshape(f(:,i),nx,ny)';
            fstar7(Top,Left) = f7(Bottom,Right);
            fstar7=reshape(fstar7',nodenums,1);
        elseif i==8
            fstar8=reshape(fstar(:,i),nx,ny)'; f8=reshape(f(:,i),nx,ny)';
            fstar8(Bottom,Left) = f8(Top,Right);
            fstar8=reshape(fstar8',nodenums,1);
        elseif i==9
            fstar9=reshape(fstar(:,i),nx,ny)'; f9=reshape(f(:,i),nx,ny)';
            fstar9(Bottom,Right) = f9(Top,Left);
            fstar9=reshape(fstar9',nodenums,1);
        end
    end
fstar = [fstar1 fstar2 fstar3 fstar4 fstar5 fstar6 fstar7 fstar8 fstar9];
end
