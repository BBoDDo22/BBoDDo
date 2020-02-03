# GuroMulgae

% Hemal Dave Project 4B: Medical Image Processing
% this is the hybrid median filter. several options will be provided to
% plot the filter on edge processing things.kernel size and pad value is
% the user input values.I am not going to document all the things because
% most of the part is repeted.
clc;
close all;
clear all;
%% Choise panel for padding method and chossing kernel size.
Kernel_size = input('Enter the kernel Size : ');% chossing the kernel size
Choise = menu('Choose Filter Technique','1) no processing of border pixels','2) padding with a fixed value','3) padding with image reflection');
if Choise == 2
    pad = input('Enter your fixed pad value (must be between 0 to 1) : ');
end
if Choise == 1
    ImageA = im2double(imread('BME7112_Data_File_4A.tif'));
    [x y z] = size(ImageA);% check if the image is RGB or grey scale if image is RGB this will convert into greyscale.
if z == 1% double check if image is rgb , will convert to grey scale
    ImageA = ImageA;
else
ImageA = rgb2gray(ImageA);
end
figure(1)
subplot(1,2,1)
imshow(ImageA)
title('Original Image');
if mod(Kernel_size,2) == 0% check the kernel size if it is even it will increment by 1
    Kernel_size = Kernel_size+1;
else
    Kernel_size = Kernel_size;
end
start = (round(Kernel_size/2));% adjustment for choosing the element of + & x.
as = Kernel_size - start+1;% adjustment for choosing the element of + & x.
As = as-1;% adjustment for choosing the element of + & x.
[x,y,z] = size(ImageA);
filtered_ImageA = ImageA;
Count = 0;
for a = 1:z
    for i = start:x-start+1
        for j = start:y-start+1
            sum = 0;
            
            for k = 1:Kernel_size
                for l = 1:Kernel_size
                    sum = sum + 1;
                    count(sum) = ImageA(k+i-as,l+j-as,a);% making the vector of the whole kernel.
                end
            end
            cross_1 = count(Kernel_size:Kernel_size-1:end-Kernel_size+1);% count the 'x' from right top corner to left bottom from the kernel
            cross_2 = count(1:Kernel_size+1:end);% count the 'x' from left top corner to right bottom from the kernel.
            cross_2(Kernel_size-As) = [];% removing the central pixel which is repeting twice while selecting pixels from x.
            cross_main = sort([cross_1 cross_2]);% sort the matrix & merge both cross matrix
            plus_1 = count(as:Kernel_size:end);% count the '+' from top to botom
            plus_2 = count((Kernel_size*As)+1:Kernel_size*(As+1));% count the '+' from left to right on central pixel line.
            plus_2(Kernel_size-As) = [];% removing the central pixel which is repeting twice while selecting pixels from +.
            plus_main = sort([plus_1 plus_2]);% sort and merge both vecters of +.
            median_cross_main = median(cross_main);% finding median for x
            median_plus_main = median(plus_main);% finding median for +.
            final_kernel = [median_cross_main,count((Kernel_size*as)-As),median_plus_main];% merge this whole median with the central pixel
            final_kernel = sort(final_kernel);% sort the vecter
            filtered_ImageA(k+i-Kernel_size,j+l-Kernel_size,a) = median(final_kernel);% find the median and replace the central pixel value with median value.
        end
    end
end
subplot(1,2,2)
imshow(filtered_ImageA)
title('Hybrid Median Filter with no processing on the edge');
figure()
subplot(1,2,2)
Acropped = filtered_ImageA(152:218,161:194);
imshow(Acropped);
title('Cropped Image after filter')
subplot(1,2,1)
Acropped1 = ImageA(152:218,161:194);
imshow(Acropped1);
title('Crooped Image before filter')

else if Choise == 2
        ImageA = im2double(imread('BME7112_Data_File_4A.tif'));
        [x y z] = size(ImageA);% check if the image is RGB or grey scale if image is RGB this will convert into greyscale.
if z == 1
    ImageA = ImageA;
else
ImageA = rgb2gray(ImageA);
end
figure(1)
subplot(1,3,1)
imshow(ImageA)
title('Original Image');
%%
if mod(Kernel_size,2) == 0
    Kernel_size = Kernel_size+1;
else
    Kernel_size = Kernel_size;
end
as = Kernel_size -1;
As = as-1;
Asd = round(Kernel_size/2)-1;
asdf = Kernel_size-Asd;
%% Padding
[x,y,z] = size(ImageA);
ImageP = ones(size(ImageA)+as);
ImageP = ImageP.*pad;% change the Pad value which is given by user.
ImageP(((as/2)+1):end-(as/2),((as/2)+1):end-(as/2)) = ImageA;
ImageA = ImageP;
subplot(1,3,2)
imshow(ImageA)
title('Padded Image');
filtered_ImageA = ImageA;
[x y z] = size(ImageA);
% Hybrid Median Filter
for a = 1:z
    for i = (as/2):x-((as/2)+1)
        for j = (as/2):y-((as/2)+1) 
            sum = 0;
            for k = 1:Kernel_size
                for l = 1:Kernel_size
                    sum = sum + 1;
                    count(sum) = ImageA(k+i-(as/2),l+j-(as/2));
                end
            end
            cross_1 = count(Kernel_size:Kernel_size-1:end-Kernel_size+1);
            cross_2 = count(1:Kernel_size+1:end);
            cross_2(Kernel_size-Asd) = [];
            cross_main = sort([cross_1 cross_2]);
            plus_1 = count((Asd+1):Kernel_size:end);
            plus_2 = count((Kernel_size*Asd)+1:Kernel_size*(Asd+1)); 
            plus_2(Kernel_size-Asd) = [];
            plus_main = sort([plus_1 plus_2]);
            median_cross_main = median(cross_main);
            median_plus_main = median(plus_main);
            final_kernel = [median_cross_main,count((Kernel_size*Asd)+asdf),median_plus_main];
            final_kernel = sort(final_kernel);
            filtered_ImageA(k+i-as,j+l-as,a) = median(final_kernel); 
        end
    end
end
filtered_ImageA =filtered_ImageA;
subplot(1,3,3)
imshow(filtered_ImageA)
title('Hybrid Median Filtered Image Padding with fixed value ')
figure()
subplot(1,2,2)
Acropped = filtered_ImageA(152:218,161:194);
imshow(Acropped);
title('Cropped Image after filter')
subplot(1,2,1)
Acropped1 = ImageA(152:218,161:194);
imshow(Acropped1);
title('Crooped Image before filter')
    else
        ImageA = im2double(imread('BME7112_Data_File_4A.tif'));
[x y z] = size(ImageA);
if z == 1
    ImageA = ImageA;
else
ImageA = rgb2gray(ImageA);
end
figure(1)
subplot(1,3,1)
imshow(ImageA)
title('Original Image');
if mod(Kernel_size,2) == 0
    Kernel_size = Kernel_size+1;
else
    Kernel_size = Kernel_size;
end
as = Kernel_size -1;
As = as-1;
Asd = round(Kernel_size/2)-1;
asdf = Kernel_size-Asd;
[x,y,z] = size(ImageA);
ImageP = ones(size(ImageA)+as);
% reflection padding.
ImageP(((as/2)+1):end-(as/2),((as/2)+1):end-(as/2)) = ImageA;
ImagePA = flipud(ImageP(((as/2)+1):(as),:));
ImageP(1:(as/2),:) = ImagePA;
ImagePB = flipud(ImageP(end-as+1:end-(as/2),:));
ImageP(end-(as/2)+1:end,:) = ImagePB;
ImagePC = fliplr(ImageP(:,((as/2)+1):(as)));
ImageP(:,1:(as/2)) = ImagePC;
ImagePD = fliplr(ImageP(:,end-as+1:end-(as/2)));
ImageP(:,end-(as/2)+1:end) = ImagePD;
filtered_ImageA = ImageP;
subplot(1,3,2)
imshow(filtered_ImageA);
title('Padded Image');
ImageA = ImageP;
[x y z] = size(ImageA);
for a = 1:z
    for i = (as/2):x-((as/2)+1)
        for j = (as/2):y-((as/2)+1) 
            sum = 0;
            for k = 1:Kernel_size
                for l = 1:Kernel_size
                    sum = sum + 1;
                    count(sum) = ImageA(k+i-(as/2),l+j-(as/2));
                end
            end
            cross_1 = count(Kernel_size:Kernel_size-1:end-Kernel_size+1);
            cross_2 = count(1:Kernel_size+1:end);
            cross_2(Kernel_size-Asd) = [];
            cross_main = sort([cross_1 cross_2]);
            plus_1 = count((Asd+1):Kernel_size:end);
            plus_2 = count((Kernel_size*Asd)+1:Kernel_size*(Asd+1)); 
            plus_2(Kernel_size-Asd) = [];
            plus_main = sort([plus_1 plus_2]);
            median_cross_main = median(cross_main);
            median_plus_main = median(plus_main);
            final_kernel = [median_cross_main,count((Kernel_size*Asd)+asdf),median_plus_main];
            final_kernel = sort(final_kernel);
            filtered_ImageA(k+i-as,j+l-as,a) = median(final_kernel);                
        end
    end
end
% Here I am doing the padding of the reflection with the updated value of
% the hybrid median filter.
ImageP = ones(size(ImageA)+as);
filtered_ImageAP = filtered_ImageA(((as/2)+1):end-(as/2),((as/2)+1):end-(as/2));
filtered_ImageA(((as/2)+1):end-(as/2),((as/2)+1):end-(as/2)) = filtered_ImageAP;
filtered_ImageAA = flipud(filtered_ImageA(((as/2)+1):(as),:));
filtered_ImageA(1:(as/2),:) = filtered_ImageAA;
filtered_ImageAB = flipud(filtered_ImageA(end-as+1:end-(as/2),:));
filtered_ImageA(end-(as/2)+1:end,:) = filtered_ImageAB;
filtered_ImageAC = fliplr(filtered_ImageA(:,((as/2)+1):(as)));
filtered_ImageA(:,1:(as/2)) = filtered_ImageAC;
filtered_ImageAD = fliplr(filtered_ImageA(:,end-as+1:end-(as/2)));
filtered_ImageA(:,end-(as/2)+1:end) = filtered_ImageAD;
subplot(1,3,3)
imshow(filtered_ImageA)
title('Hyber median filtered Image padding with reflection');
figure()
subplot(1,2,2)
Acropped = filtered_ImageA(152:218,161:194);
imshow(Acropped);
title('Cropped Image after filter')
subplot(1,2,1)
Acropped1 = ImageA(152:218,161:194);
imshow(Acropped1);
title('Crooped Image before filter')
    end
end


