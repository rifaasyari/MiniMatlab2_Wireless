function [ output_args ] = time_freq_analysis( S, threshold)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

S_binary = S > threshold ; 

%%
clc; clear all ; close all; 
%dilation neighborhood
dil = strel('square',2 )



er = strel('square',2)

%erosion neighborhood 


%
A =   [ zeros(1,12 ) ; zeros(1,5 ), 1 ,zeros(1,5) 1 ;  ... 
        0 , 1 ,zeros(1,2) , ones(1,1) , zeros(1,4) , 1 zeros(1,2) ; ...
        ones(1,3) , zeros(1,5), ones(1,4)   ; ... 
        ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
        ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
        zeros(1,12) ; ...
        zeros(1,5) , 1, zeros( 1,4), 1 , 0     ;   ... 
        zeros(1,12) ; ...
        zeros(1,12) ; ...
        zeros(1,12) ; ...
        ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
        1, 0 , 1 ,  zeros(1,5), 1, 0 , 1, 0 ; ...
        ones(1,3) , zeros(1,5), ones(1,3), 0 ; ...
        zeros(1,12) ; ...
        zeros(1,12) ; ...
        zeros(1,12) ; ...
        zeros(1,12) ; ...
        
        zeros(1,8), ones(1,3), 0 ; ...
        zeros(1,4),1, zeros(1,3), 1 , 1, 0 , 0 ; ...
        zeros(1,8), 0 , 1, 0 , 0 ;...   
        zeros(1,12) ; ...
        zeros(1,12) ; ...
         ] ;

%
 A = padarray(A, [ 4,4]);
%I2 = imdilate(I1,dil) ;
%
I1 = imdilate(A,dil)  ;
I2 = imerode(I1,er) ;
%


% MEDIAN FILTER REMOVES SALT AND PEPPER NOISE
B = medfilt2(I2, [3,3 ]);




figure; imshow(A,'InitialMagnification', 800); 
figure; imshow(I2,'InitialMagnification', 800); 
figure; imshow(B,'InitialMagnification', 800); 

D = bwlabel(B);

measurements = regionprops(D,'Centroid');

measurements.Centroid

horizontal_projection = sum(B , 1) ;

max_h = max(horizontal_projection) ; 
horizontal_coords = find( horizontal_projection >= 0.1*max_h) ; 

vertical_projection = sum(B, 2) ; 
%%
 figure; plot(horizontal_projection) ;
  figure; plot(vertical_projection) ;

end

