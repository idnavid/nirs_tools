function [raw_data,hb_concentrations,source_detector_labels,fs,brain_map] = load_data()
addpath(genpath('C:\Users\nshokouhi\Downloads\NIRS\HRF_estimation\debug_code\nirs_spm_dir\spm12'))
addpath(genpath('../tools/'))
% load raw NIRS data. These contain optical densities.
fs = 12.5;
source_detector_labels = active_channels();

brain_map = label2map(source_detector_labels,3,5);
filelist = '../data/data_list.txt';
fin = fopen(filelist);
C = textscan(fin,'%s');
raw_data = cell(length(C{1}),2);
hb_concentrations = cell(length(C{1}),3);% hbO,hbR,hbT
for i = 1:length(C{1})
    dirname = C{1}{i};
    dirname = ['../',dirname];
    split1 = strsplit(dirname,'/');
    basename = split1{end};
    wl1_filename = [dirname,'/NIRS-',basename,'.wl1']
    wl2_filename = [dirname,'/NIRS-',basename,'.wl2']
    raw_data(i,1) = read_file(wl1_filename);
    raw_data(i,2) = read_file(wl2_filename);
    
    % calculate hb responses
    filename = [dirname,'/NIRS-',basename];
    [y, P] = read_nirscout(filename);
    [X, P] = spm_fnirs_calc_od(y, P);
    [hbo, hbr, hbt] = spm_fnirs_calc_hb(X, P);
    hb_concentrations(i,1)={hbo};
    hb_concentrations(i,2)={hbr};
    hb_concentrations(i,3)={hbt};
end
end


%% Label active data channels
function [labels] = active_channels()
% The channel mask below shows which source-detector
% pairs were activated. It was obtained from the
% NIRS header file (*.hdr).
% rows correspond to sources and columns to detectors. 
channel_mask = [0       0       1       1       1;
                1       0       1       1       0;
                0       1       1       0       1];
labels = logical(channel_mask(:));
end

%% Find source-detector location
function brain_map = label2map(labels,num_sources,num_detectors)
% converts channel-mask to coordinates of source-detector pairs
% 
% The main problem with channel indexes is that the physical location 
% of channels is unclear. label2map takes the labels and maps the desired
% channels (those with 1 in their corresponding entry) to a physical
% location in the MNI coordinate. 
% Inputs:
%       labels:     mask of 1s and 0s for source-detector pairs.
%                   sources are rows, detectors are columns
%
% Output:
%       brain_map: cell of pairs containing source and detector locations.
%                  in each cell the first row is the source and the second
%                  row is the detector. Each cell corresponds to a channel.
                    
labels = double(reshape(labels,num_sources,num_detectors));
brain_map = cell(sum(labels(:)),1);

% NOTE: locations are hard-coded and are specific to our montage. 
%       The file ".txt" was used to find the x,y,z coordinates. 

% Locations: x, y , z
source_locs = [0	0	88;
               -52	-52	47;
               -52	52	47];

detector_locs = [-83	-27 -3;
                 -83	27	-3;
                 -63	0	61;
                 0      -63	61;
                 0      63	61];


channel_idx = 1;
for detector_idx = 1:size(labels,2)
    for source_idx = 1:size(labels,1)
        if labels(source_idx,detector_idx)==1
            brain_map{channel_idx} = [source_locs(source_idx,:);
                detector_locs(detector_idx,:)];
            channel_idx = channel_idx + 1;
        end
    end
end
end

%% Remove NaNs from data file
function data = read_file(filename)
% reads text file. 
temp = dlmread(filename);
temp(isnan(temp)) = 0;
data = {temp};
end