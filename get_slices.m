function [fixed, R_fixed, moving, R_moving, roi, R_roi] = get_slices(fixed_file, moving_file, roi_file, reg_dimensions, modality_tag, visit_nb)

% A script to load and pre-process the fixed image, moving image and ROI to register.
% Input: fixed/moving/roi image complete filepaths (string), dimensions of registration (2 or 3), modality number (1 to 5, see main.m)
% Output: correctly formatted images and corresponding spatial references

%% Prepare parameters depending on registration type and data type

% recall indices are as follows:
% raw_options = {'T2* STARMAP', 'DWI', 'MOLLI', 'VFA', 'DCE'};
% map_options = {'T2* map', 'ADC map', 'T1 map (MOLLI)', 'T1 map (VFA)', 'k-trans map'};

% 4th dimension timepoint for fixed image (index 1 for T2* STARMAP, 2 for DWI, 3 for MOLLI etc)
% 3rd dimension slice selection for fixed image
timepoint_4d = [0, 1, 9/11, 0, 2/65];
timepoint_3d = [0, 8, 0, 8, 14];  

% slice number for the T2-weighted image (moving image)
% **** TAKE CARE, DIFFERENT SLICE FOR EACH VISIT ****
if visit_nb == 1
    slice_moving = 5;
elseif visit_nb==2
    slice_moving = 10;
end

%% Process data to obtain correct 2D or 3D format

% "fixed" is the selected raw data reference map
fixed = niftiread(fixed_file);
fixed = fliplr(rot90(fixed,3));
fixed_info = niftiinfo(fixed_file);
% here, select correct timepoint if data is 4d
if length(size(fixed)) == 4
    t = round( timepoint_4d(modality_tag) * size(fixed,4) );
    fixed = fixed(:,:,:,t);    % last index is timepoint
end
% add padding if necessary
if size(fixed, 3)~= 1 && size(fixed, 3) < 16
    padding = round((16 - size(fixed,3))/2);
    padded_fixed = uint16(zeros(256, 256, size(fixed,3)+2*padding));
    padded_fixed(:, :, padding+1:size(padded_fixed,3)-padding) = fixed;
    fixed = padded_fixed;
end

% "moving" is the selected input image (T2 weighted)
moving = niftiread(moving_file);
moving = fliplr(rot90(moving,3));
moving_info = niftiinfo(moving_file);
% add padding if necessary
if size(moving, 3)~= 1 && size(moving, 3) < 16
    padding = round((16 - size(moving,3))/2);
    padded_moving = uint16(zeros(256, 256, size(moving,3)+2*padding));
    padded_moving(:, :, padding+1:size(padded_moving,3)-padding) = moving;
    moving = padded_moving;
end

% "roi" is the region of interest in the input image (T2-weighted)
roi = niftiread(roi_file);
roi = fliplr(rot90(roi,3));
roi_info = niftiinfo(roi_file);
% add padding if necessary
if size(roi, 3)~= 1 && size(roi, 3) < 16
    padding = round((16 - size(roi,3))/2);
    padded_roi = uint16(zeros(256, 256, size(roi,3)+2*padding));
    padded_roi(:, :, padding+1:size(padded_roi,3)-padding) = roi;
    roi = padded_roi;
end

% get spatial references using the get_relative_imrefs.m script
[R_fixed, R_moving] = get_relative_imrefs(fixed_info, moving_info, reg_dimensions);
[~, R_roi] = get_relative_imrefs(fixed_info, roi_info, reg_dimensions);

% extract individual slices for 2D registration
%  /!\ manual for now but could be automated /!\
if reg_dimensions == 2
    if length(size(fixed)) == 3
        slice_fixed = timepoint_3d(modality_tag);
        fixed = fixed(:,:,slice_fixed);
    end
    moving = moving(:,:,slice_moving);
    roi = roi(:,:,slice_moving);
end

% ---- END OF PREPROCESSING ----

end 

