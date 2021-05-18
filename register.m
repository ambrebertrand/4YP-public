function reg_roi = register(fixed, R_fixed, moving, R_moving, roi, R_roi, modality_tag, visit_nb, display_flag)

% A script to load the fixed image, moving image and ROI to register.
% Input: fixed/moving/roi image complete filepaths (string), dimensions of registration (2 or 3), modality number (1 to 5, see main.m)
% Output: registered ROI as a 2d binary mask

%% Performs a 2D-to-2D image registration of structural to functional raw data
% Or 3D to 3D if both images are 3D

% Select optimiser, metric and registration type
[optimizer,metric] = imregconfig('multimodal');

% Parameters to be tuned for specific application
optimizer.InitialRadius = 0.003;         % 0.01
optimizer.Epsilon = 1.5e-4;               % 5e-6
optimizer.GrowthFactor = 1.01;         % 1.01
optimizer.MaximumIterations = 5000;         % 5000

% Compute transform matrix
tform = imregtform(moving, R_moving, fixed, R_fixed, 'affine', optimizer, metric);

% ROI is transformed using transformation matrix "tform" computed above
[reg_moving, R_reg_moving] = imwarp(moving, R_moving, tform, 'OutputView', R_fixed);

% use R_moving here again as the ROI spatial refs can be dodgy 
[reg_roi, R_reg_roi] = imwarp(roi,R_moving,tform,'OutputView',R_fixed);


% ----------- Visualisations -----------

% recall indices are as follows:
raw_options = {'T2* STARMAP', 'DWI', 'MOLLI', 'VFA', 'DCE'};
map_options = {'T2* map', 'ADC map', 'T1 map (MOLLI)', 'T1 map (VFA)', 'k-trans map'};

% displays
rawmod_name = cell2mat(raw_options(modality_tag));
ArgsForDisplay = {fixed, R_fixed, moving, R_moving, reg_moving, R_reg_moving, roi, R_roi, reg_roi, R_reg_roi, visit_nb, rawmod_name};
if display_flag == 'y'
    display_results(ArgsForDisplay);
end 

% %Show fixed and moving images overlaid before and after registration
% figure
% subplot(1,2,1)
% imshowpair(fixed, R_fixed, moving, R_moving, 'Scaling','independent')
% title('Before registration')
% axis off
% subplot(1,2,2)
% imshowpair(fixed, R_fixed, reg_moving, R_reg_moving,'Scaling','independent');
% title('After registration')
% hold on
% axis off
% 
% % Show original and registered ROI
% 
% figure
% h1 = axes;
% imagesc(moving);
% colormap(h1,gray)
% set(h1,'ydir','reverse');
% hold on
% axis equal
% axis off
% title('ROI on original moving image')
% h2 = axes;
% imagesc('CData', moving.*roi, 'AlphaData', roi, 'AlphaDataMapping', 'scaled');
% set(h2,'color','none','visible','off')
% colormap(h2,parula)
% set(h2,'ydir','reverse');
% %   caxis([0 20])
% c=colorbar;
% c.Location = 'east';
% c.AxisLocation = 'out';
% linkaxes([h1 h2])
% set(h1, 'xlim', [0 size(moving, 2)], 'ylim', [0 size(moving, 1)])
% axis equal
% 
% figure
% h1 = axes;
% imagesc(fixed);
% colormap(h1,gray)
% set(h1,'ydir','reverse');
% hold on
% axis equal
% axis off
% title('Registered ROI on fixed image')
% h2 = axes;
% imagesc('CData', int16(fixed).*reg_roi, 'AlphaData', reg_roi, 'AlphaDataMapping', 'scaled');
% set(h2,'color','none','visible','off')
% colormap(h2, parula)
% set(h2,'ydir','reverse');
% %   caxis([0 20])
% c=colorbar;
% c.Location = 'east';
% c.AxisLocation = 'out';
% linkaxes([h1 h2])
% set(h1, 'xlim', [0 size(fixed, 2)], 'ylim', [0 size(fixed, 1)])
% axis equal

% ----------- end of visualisation -----------

end

