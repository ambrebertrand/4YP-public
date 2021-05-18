%% Main script for REACT app (Registration, Extraction and Analysis of Cancer Tumour MRI data) 
% Ambre Bertrand, University of Oxford, May 2021

% ---- Read carefully before running the scripts! ----

% 1:  MATLAB scipt dependencies (make sure everything is in the same CODE_PATH):
% get_relative_imrefs.m
% get_slices.m
% register.m
% display_results.m 
% p_call.m
% histograms.m

% 2:  Make sure to change paths set by default to Ambre's computer. The 
% correct code is commented out, uncomment that and remove Ambre's paths.

% 3:  If dealing with non-ART trial data, correct volume slices might
% need to be checked again. These are set at the start of the get_sices.m
% script in array form for the functional modalities and a single number
% for the T2-weighted image which changes from one visit to the next. 

% ------------

clear
close all

%% Directories (default and user inputs)

% ---- These are the default paths (Ambre's computer) ----
CODE_PATH = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Code/my_matlab';
p0 = genpath(CODE_PATH);
DATA_PATH_1 = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI250717/';
p1 = genpath(DATA_PATH_1);
DATA_PATH_2 = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI180817/';
p2 = genpath(DATA_PATH_2);
% -------------------------------------------------------------------

%% Dialog box for user to change directory paths if required
prompt = {'Enter path with .m files:', 'Enter directory path for Visit 1:', 'Enter directory path for Visit 2:'};
default = {CODE_PATH, DATA_PATH_1, DATA_PATH_2};
dlgtitle = 'Input';
dims = [1 100];
answer = inputdlg(prompt,dlgtitle,dims, default);
CODE_PATH = cell2mat(answer(1));
DATA_PATH_1 = cell2mat(answer(2));
DATA_PATH_2 = cell2mat(answer(3));


%% Set up imaging modality options

% Indices in cell correspond to indices in checkboxes
raw_options = {'T2* STARMAP', 'DWI', 'MOLLI', 'VFA', 'DCE'};
map_options = {'T2* map', 'ADC map', 'T1 map (MOLLI)', 'T1 map (VFA)', 'k-trans map'};

% Checkbox figure to select modalities of interest
h.f = figure('units','pixels','position',[300,300,200,250],...
             'toolbar','none','menu','none');
annotation('textbox', [0,0.99,1,0], 'string', 'Select imaging modalities for input.')

% create checkboxes for the modalities
for i=1:length(raw_options)
    label = raw_options{i};
    h.c(i) = uicontrol('style','checkbox','units','pixels',...
                'position',[10,(225-30*i),150,15],'string',label);
end
            
% create button to  confirm selection  
h.p = uicontrol('style','pushbutton','units','pixels',...
                'position',[40,5,70,20],'string','OK');

% callback returns array 'checked' with modality numbers selected
% e.g. checked = [1 3 4] if you've picked t2*, molli, oe molli
set(h.p, 'callback', @(src, event) p_call(src, event, h));
waitfor(h.f);
nb_of_modalities = size(checked,2);

%% Initialise cells to store variables

nb_of_options = length(raw_options);

% Create cells to store outputs
roi_vals_1 = cell(1,nb_of_options);
roi_vals_2 = cell(1,nb_of_options);
roi_stats_1 = cell(1,nb_of_options);
roi_stats_2 = cell(1,nb_of_options);
reg_roi_1 = cell(1,nb_of_modalities);
reg_roi_2 = cell(1,nb_of_modalities);

% select registration type (*** ONLY WORKS FOR 2 as of 01/02/21 ***)
prompt = {'Enter registration dimensions (2 or 3):'};
dlgtitle = 'Input';
dims = [1 100];
answer = inputdlg(prompt,dlgtitle,dims);
reg_dimensions = str2double(cell2mat(answer(1)));


%% Perform registration and ROI analysis for each modality

% loop for visit 1 and 2
for v=1:2
    % change to directory of correct visit data
    if v==1
        addpath(p1);
        cd(DATA_PATH_1)
    else
        addpath(p2);
        cd(DATA_PATH_2)
    end
    
    % select T2 weighted image and ROI for that visit (same for all fixed)
    f = msgbox(sprintf('Select T2-weighted image (moving) for Visit %d', v)); uiwait(f);
    moving_file = uigetfile('*.nii.gz','Select moving image');
    f = msgbox(sprintf('Select ROI for Visit %d', v)); uiwait(f);
    roi_file = uigetfile('*.nii.gz','Select ROI');
    
%     % ---- use this for quick testing ----
%     if v==1
%         moving_file = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI250717/6_ax_frfse_t2_hr.nii.gz';
%         roi_file = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI250717/0_RTSTRUCT_MRI_1_1-label.nii.gz';
%         roi_file_1 = roi_file;
%     else 
%         moving_file = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI180817/4_ax_frfse_t2_hr.nii.gz';
%         roi_file = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/Data/ART_128_MRI180817/0 RTSTRUCT MRI_1_2-label.nii.gz';
%         roi_file_2 = roi_file;
%     end    
%     % --------------------
%     
    % create a list to store file name of each required fixed modality
    modality_idx = zeros(1,nb_of_modalities);
    fixed_files = cell(1,nb_of_modalities);
    
    % load all selected RAW fixed modalities to register
    for m=1:nb_of_modalities
        modality_idx(m) = checked(m);
        modality_name = cell2mat(raw_options(modality_idx(m)));
        f = msgbox(sprintf('Select %s (fixed) for Visit %d', modality_name, v)); uiwait(f);
        fixed_files{m} = uigetfile('*.nii.gz');        
    end 
   
    % preprocessing and registration for each fixed modaity
    for m=1:nb_of_modalities               
        addpath(p0);
        cd(CODE_PATH)   
           
        % load and process images selected to get correct 2D slices
        [fixed, R_fixed, moving, R_moving, roi, R_roi] = get_slices(fixed_files{m}, moving_file, roi_file, reg_dimensions, modality_idx(m), v);
        
        % include registration visualisations? ('y' if yes, otherwise no)
        display_flag = 'no';
        % compute registered ROI binary mask
        reg_roi = register(fixed, R_fixed, moving, R_moving, roi, R_roi, modality_idx(m), v, display_flag);
        reg_roi = fliplr(rot90(reg_roi,2));
        
        % assign computed ROI to correct modality
        if v==1
            reg_roi_1{m} = reg_roi;
        else
            reg_roi_2{m} = reg_roi;
        end  
    end
    
  
    %% Apply binary ROI to quantitative data map to extract ROI values
    
    % reset lists to store indices and file names
    quantimap_file = cell(1,nb_of_modalities);
     
    % apply roi to each quantitative map as selected
    for m=1:nb_of_modalities

        if v==1
            addpath(p1);
            cd(DATA_PATH_1)
        else
            addpath(p2);
            cd(DATA_PATH_2)
        end

        % select quantitative map to apply ROI mask to
        modality_name = cell2mat(map_options(modality_idx(m)));
        f = msgbox(sprintf('*** Quantitative map input *** \n Select %s for Visit %d', modality_name, v)); uiwait(f);
        quantimap_file{m} = uigetfile('*.nii.gz');

        % --------- similar preprocessing as in register.m ----------------
        quanti = niftiread(quantimap_file{m});
        % rotate the non-T1 maps correctly
        if (modality_idx(m) == 1) || (modality_idx(m) == 2) || (modality_idx(m) == 5)   
            quanti = rot90(quanti);              
        else
            quanti = fliplr(rot90(quanti,3));    
        end
        quanti_info = niftiinfo(quantimap_file{m});
        % select correct timepoint if data is 4d
        timepoints = [13/16, 1, 6/11, 0, 2/65];
        if length(size(quanti)) == 4
            if modality_idx(m)~=4     % for everyone but VFA
                t = round(timepoints(modality_idx(m)) * size(quanti,4));
            end
            quanti = quanti(:,:,:,t);    % last index is timepoint
        end
        % add padding if necessary
        if size(quanti, 3)~= 1 && size(quanti, 3) < 16
            padding = round((16 - size(quanti,3))/2);
            padded_quanti = uint16(zeros(256, 256, size(quanti,3)+2*padding));
            padded_quanti(:, :, padding+1:size(padded_quanti,3)-padding) = quanti;
            quanti = padded_quanti;
        end  
        % extract individual slice for 2D case
        slice_nb=10;     % same as for registration
        if reg_dimensions == 2
            if length(size(quanti)) == 3
                quanti = quanti(:,:,slice_nb);
            end
        end
        % ----------------------------------------------------

        % extract tumour values using ROI
        if v==1
            quantimap1 = int16(quanti);
            roi_vals_1{modality_idx(m)} = reg_roi_1{m} .* quantimap1;
        else
            quantimap2 = int16(quanti);
            roi_vals_2{modality_idx(m)} = reg_roi_2{m} .* quantimap2;
        end

        % visualisations
        if v==1
            figure
            h1 = axes;
            if modality_idx(m)==1
                imagesc(quantimap1,[0,300])
            else
                imagesc(quantimap1);
            end
            colormap(h1,gray)
            set(h1,'ydir','reverse');
            hold on
            axis equal
            axis off
            title(sprintf('ROI on %s, Visit 1', modality_name))
            h2 = axes;
            imagesc('CData', quantimap1.*reg_roi_1{m}, 'AlphaData', reg_roi_1{m}, 'AlphaDataMapping', 'scaled');
            set(h2,'color','none','visible','off')
            colormap(h2,parula)
            set(h2,'ydir','reverse');
            %   caxis([0 20])
            c=colorbar;
            c.Location = 'east';
            c.AxisLocation = 'out';
            linkaxes([h1 h2])
            set(h1, 'xlim', [0 size(quantimap1, 2)], 'ylim', [0 size(quantimap1, 1)])
            axis equal
        else
            figure
            h1 = axes;
            if modality_idx(m)==1
                imagesc(quantimap2,[0,300])
            else
                imagesc(quantimap2);
            end
            colormap(h1,gray)
            set(h1,'ydir','reverse');
            hold on
            axis equal
            axis off
            title(sprintf('ROI on %s, Visit 2', modality_name))
            h2 = axes;
            imagesc('CData', quantimap2.*reg_roi_2{m}, 'AlphaData', reg_roi_2{m}, 'AlphaDataMapping', 'scaled');
            set(h2,'color','none','visible','off')
            colormap(h2, parula)
            set(h2,'ydir','reverse');
            %   caxis([0 20])
            c=colorbar;
            c.Location = 'east';
            c.AxisLocation = 'out';
            linkaxes([h1 h2])
            set(h1, 'xlim', [0 size(quantimap2, 2)], 'ylim', [0 size(quantimap2, 1)])
            axis equal  
        end
    end

end

  
      
%% ROI statistics and histograms

% create array to store volume for visits 1 and 2
tumour_volume = zeros(1,2);

% TUMOUR VOLUME
for v=1:2
    if v==1
        roi_file = roi_file_1;
    else 
        roi_file = roi_file_2;
    end
    roi_data = niftiread(roi_file);
    roi_info = niftiinfo(roi_file);
    dims = roi_info.PixelDimensions;
    tumour_voxels = sum(roi_data(:));
    % compute volume in cm3
    tumour_volume(v) = tumour_voxels * dims(1) * dims(2) * dims(3) * 0.001;
end
    

%% HISTOGRAMS - loop stat analysis over number of modalities selected
units = {'T2* (ms)', 'ADC (x10^{-6} mm^2/s)', 'T1 (ms)', 'T1 (ms)', 'k^{trans} (min^{-1})'};
% initialise cell to store stats for each modality
csv_outputs = cell(1,nb_of_modalities);
for m=1:nb_of_modalities

    % select correct modality for each loop
    modality_idx = checked(m);
    modality_name = cell2mat(map_options(modality_idx));
    modality_units = cell2mat(units(modality_idx));
    
    % output: numvoxi, meani, stdi, mixingproportionsi, mui, sigmai (i=1,2)
    % ... change in volume, change in mean, change in lowpeak
    [output,values_for_csv] = histograms(roi_vals_1{modality_idx}, roi_vals_2{modality_idx}, modality_name, modality_units);
    csv_outputs{m} = values_for_csv;
    
    % get individual values for visits 1 and 2
    blob = num2cell(output(1:3));   % just a dummy for the next line
    [numvox1, mean1, std1] = deal(blob{:});
    prop1 = output(4:6);
    mu1 = output(7:9);
    sigma1 = output(10:12);

    blob = num2cell(output(13:15));
    [numvox2, mean2, std2] = deal(blob{:});
    prop2 = output(16:18);
    mu2 = output(19:21);
    sigma2 = output(22:24);

    deltavol = output(25);
    deltamean = output(26);
    deltalowpeak = output(27);
    deltamu = 100*(mu2-mu1)./mu1;

    % below, record the values that we want to keep for final outputs
    roi_stats_1{modality_idx} = [numvox1, mean1, std1];
    roi_stats_2{modality_idx} = [numvox2, mean2, std2];

%     % plot bar graphs
%     figure
%     labels = {'Visit 1 GMM proportions (%)';'Visit 2 GMM proportions (%)';'Change in mean (%)'};
%     y = [100*prop1(:)'; 100*prop2(:)';deltamu(:)'];
%     b = bar(y);
%     set(gca,'xticklabel',labels)
%     title('Tumour values GMM components')
%     ylim([-100 100])
% 
%     for i=1:3
%         xtips = b(i).XData;
%         ytips = b(i).YData;
%         labels = string(b(i).YData);
%         if i==1
%             align='right';
%         elseif i==2
%             align='center';
%         else
%             align='left';
%         end
%         text(xtips,ytips,labels,'HorizontalAlignment',align,'VerticalAlignment','bottom')
%     end
end


%% Basic statistical plots

% input potential healthy/diseased threshold 
thresholds = [1000,1500,1150,1000,1000];

% BAR CHARTS
for m=1:nb_of_modalities
    
    % again, select correct modality for each loop
    modality_idx = checked(m);
    modality_name = cell2mat(map_options(modality_idx));
    
    % get mean and std values
    mean1 = roi_stats_1{modality_idx}(2);
    mean2 = roi_stats_2{modality_idx}(2);
    std1 = roi_stats_1{modality_idx}(3);
    std2 = roi_stats_2{modality_idx}(3);
    
    % bar chart with error bars
    subplot(1, nb_of_modalities, m);
    labels = {'Visit 1';'Visit 2'};
    x = 1:2;
    y = [mean1 mean2]';
    errors = [std1 std2];
    b = bar(y,0.5);
    hold on
    er = errorbar(x,y,errors,errors);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    set(gca,'xticklabel',labels);
    title(sprintf('Mean and std \n for %s', modality_name));
    ylim([0 1.2*max([mean1,mean2])]);   
    hold on 
    yline(thresholds(modality_idx), '-r', 'LineWidth',3);
end

%% Create csv with relevant stats

csv_titles = [{'Visit'},'Tumour_volume_cm3','Biomarker','Mean','SD','Number_of_GMM_components','Comp_1_mean','Comp_2_mean','Comp_3_mean','Comp_1_prop','Comp_2_prop','Comp_3_prop'];
csv_data = cell(nb_of_modalities*2, length(csv_titles));

for m=1:nb_of_modalities
    
    % select correct modality at each iteration
    modality_idx = checked(m);
    modality_name = cell2mat(map_options(modality_idx));
    
    % get data prepared from histograms script and fill in the gaps
    mod_data = num2cell(csv_outputs{1});
    mod_data{1,2} = tumour_volume(1);
    mod_data{1,3} = modality_name;
    mod_data{2,2} = tumour_volume(2);
    mod_data{2,3} = modality_name;
    
    % assign modality data to correct rows in the csv table
    csv_data((2*m-1):(2*m),:) = mod_data;  
end 

% ---uncomment below to save data to csv---

% % create cell array with column titles and data
% csv_array = vertcat(csv_titles, csv_data);
% 
% % Convert cell to a table and use first row as variable names
% csv_table = cell2table(csv_array(2:end,:),'VariableNames',csv_array(1,:));
%  
% % Write the table to a CSV file
% output_filename = '/Users/ambrebertrand/Documents/Oxford_Engineering/4th_year/4YP/ART129_data.csv';
% writetable(csv_table, output_filename)
% 






