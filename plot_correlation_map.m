function plot_correlation_map(correlations,brain_map)
% Plots the correlations between a reference channel and all the other 
% channels on a brain. 
% 
% Inputs:
%       correlation: array of correlation coefficients calculated between
%                    the signal from a reference channel and signals from
%                    all other channels. 
%       brain_map:   cell of source-detecor positions (x,y,z) coordinates
%                    for each channel. 

% NOTE: the baseline image is 2D, but brain_map has z coordinates as well,
% so it can potentially be extended to 3D brain images. 
background_image = double(rgb2gray(imread('head.png')));

% NOTE: offset was found based on trial an error and depends on the
% background image. 
offset = 116;


X = 0.2*(background_image>0); 
for i = 1:length(brain_map)
    source_pos = offset+brain_map{i}(1,:);
    detector_pos = offset+brain_map{i}(2,:);
    X(min(source_pos(1),detector_pos(1)):max(source_pos(1),detector_pos(1)),...
      min(source_pos(2),detector_pos(2)):max(source_pos(2),detector_pos(2))) = correlations(i);
end

X_filtered = 10*imgaussfilt(X,8);
X_filtered = X_filtered(1:size(background_image,1),1:size(background_image,2));
X_filtered = X_filtered.*(background_image>0);
imagesc(imrotate(X_filtered,90))
hold on 
for i = 1:length(brain_map)
    scatter(offset+brain_map{i}(2,1),offset+brain_map{i}(2,2),'b')
    scatter(offset+brain_map{i}(1,1),offset+brain_map{i}(1,2),'r')
end
