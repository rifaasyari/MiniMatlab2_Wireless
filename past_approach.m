function [PSD_binary,B ]  = past_approach( S, threshold, fs, str)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%t = 0:0.01:10;
%fs = 100; 
%x = sin(2*pi* 10* t) + 0.5 * sin(2*pi*5*t) ; 
[s f t p] = spectrogram(S,fs);
PSD = 10*log10(p);
f = f*fs/pi ; 
f1 = figure;
surf(t,f,PSD,'edgecolor','none'); axis tight; 
view(0,90); 
colorbar; 
xlabel('Time (s) ','Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
zlabel('PSD (t,f)','Fontsize',12);

title(['Original PSD', ' ' ,'(',str,')'],'Fontsize',16);  


saveas(f1,strcat(str, '_Original_PSD'),'png');


%%
%threshold = 0.1*max(PSD(:)) ; 
PSD_binary = double(PSD > threshold);


%%
%clc; clear all ; close all; 
%dilation neighborhood
dil = strel('square',2 )



er = strel('square',2)

%erosion neighborhood 


%
% A =   [ zeros(1,12 ) ; zeros(1,5 ), 1 ,zeros(1,5) 1 ;  ... 
%         0 , 1 ,zeros(1,2) , ones(1,1) , zeros(1,4) , 1 zeros(1,2) ; ...
%         ones(1,3) , zeros(1,5), ones(1,4)   ; ... 
%         ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
%         ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
%         zeros(1,12) ; ...
%         zeros(1,5) , 1, zeros( 1,4), 1 , 0     ;   ... 
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
%         ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
%         1, 0 , 1 ,  zeros(1,5), 1, 0 , 1, 0 ; ...
%         ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
        
%         zeros(1,8), ones(1,3), 0 ; ...
%         zeros(1,4),1, zeros(1,3), 1 , 1, 0 , 0 ; ...
%         zeros(1,8), 0 , 1, 0 , 0 ;...   
%         zeros(1,12) ; ...
%         zeros(1,12) ; ...
%          ] ;

% %
%  A = padarray(A, [ 4,4]);
% %I2 = imdilate(I1,dil) ;
% %
% I1 = imdilate(A,dil)  ;
% I2 = imerode(I1,er) ;
% %


I1 = imdilate(PSD_binary,dil)  ;
B = imerode(I1,er);





%figure; imshow(A,'InitialMagnification', 800); 
f2 = figure; imagesc(t,f,PSD_binary); 
title(['PSD Threshold Binary Image', ' ', 'T =','',num2str(threshold),' ' ,'(',str,')'],'Fontsize',12);
ylabel('Frequency, (Hz)','Fontsize',12);
xlabel('Time (s)','Fontsize',12); 
colorbar;
saveas(f2,strcat(str, 'PSD_Threshold_Binary_Image', '_', 'T =','_',num2str(threshold)),'png');

f3 = figure; imagesc(t,f,B);  
title(['PSD Dilation and Erosion', ' ', 'T =','',num2str(threshold), ' ' ,'(',str,')'],'Fontsize',12);
ylabel('Frequency, (Hz)','Fontsize',12);
xlabel('Time (s)','Fontsize',12);
colorbar;
saveas(f3,strcat(str,'PSD_Median_Filtered_Binary_Image', '_', 'T =','_',num2str(threshold)),'png');

D = bwlabel(B);

measurements = regionprops(D,'Centroid');

centds = measurements.Centroid ;

horizontal_projection = sum(B , 2) ;

max_h = max(horizontal_projection) ; 
horizontal_coords = find( horizontal_projection >= 0.1*max_h) ; 



end

