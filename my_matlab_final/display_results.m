function display_results(ArgsForDisplay)

[fixed, R_fixed, moving, R_moving, reg_moving, R_reg_moving, roi, R_roi, reg_roi, R_reg_roi, visit_nb, rawmod_name] = ArgsForDisplay{:};

% 2D registration results
if size(fixed, 3) == 1
    % Change the contrast of images if necessary
    % fixed = imadjust(fixed, stretchlim(fixed), []);
    % moving = imadjust(moving,stretchlim(moving), []);
    % reg_moving = imadjust(reg_moving, stretchlim(reg_moving, []);

    %Show fixed and moving images overlaid before and after registration
    figure
    subplot(1,2,1)
    imshowpair(fixed, R_fixed, moving, R_moving, 'Scaling','independent')
    title(sprintf('Before registration to %s V%d',rawmod_name, visit_nb))
    axis off
    subplot(1,2,2)
    imshowpair(fixed, R_fixed, reg_moving, R_reg_moving,'Scaling','independent');
    title(sprintf('After registration to %s V%d',rawmod_name, visit_nb))
    hold on
    axis off
    
    % rotate to get the correct view when looking at the raw images and ROIs
    fixed = fliplr(rot90(fixed,2));
    moving = fliplr(rot90(moving,2));
    roi = fliplr(rot90(roi,2));
    reg_roi = fliplr(rot90(reg_roi,2));
    
    % Show original and registered ROI       
    figure
    h1 = axes;
    imagesc(moving);
    colormap(h1,gray)
    set(h1,'ydir','reverse');
    hold on
    axis equal
    axis off
    title(sprintf('ROI on original T2 weighted image V%d',visit_nb))
    h2 = axes;
    imagesc('CData', moving.*roi, 'AlphaData', roi, 'AlphaDataMapping', 'scaled');
    set(h2,'color','none','visible','off')
    colormap(h2,parula)
    set(h2,'ydir','reverse');
%   caxis([0 20])
    c=colorbar;
    c.Location = 'east';
    c.AxisLocation = 'out';
    linkaxes([h1 h2])
    set(h1, 'xlim', [0 size(moving, 2)], 'ylim', [0 size(moving, 1)])
    axis equal
    
    figure
    h1 = axes;
    imagesc(fixed);
    colormap(h1,gray)
    set(h1,'ydir','reverse');
    hold on
    axis equal
    axis off
    title(sprintf('Registered ROI on fixed %s V%d', rawmod_name, visit_nb))
    h2 = axes;
    imagesc('CData', int16(fixed).*reg_roi, 'AlphaData', reg_roi, 'AlphaDataMapping', 'scaled');
    set(h2,'color','none','visible','off')
    colormap(h2, parula)
    set(h2,'ydir','reverse');
%   caxis([0 20])
    c=colorbar;
    c.Location = 'east';
    c.AxisLocation = 'out';
    linkaxes([h1 h2])
    set(h1, 'xlim', [0 size(fixed, 2)], 'ylim', [0 size(fixed, 1)])
    axis equal

% 3D registration results
elseif size(fixed, 3) ~= 1
    nbslice=size(fixed, 3);
    moving_original = imwarp(moving, R_moving, affine3d, 'OutputView', R_fixed);
    for i=7:10
        fixed_slice = imadjust(fixed(:,:,i), stretchlim(fixed(:,:,i)), []);
        moving_original_slice = imadjust(moving_original(:,:,i), stretchlim(moving_original(:,:,i)), []);
        reg_moving_slice = imadjust(reg_moving(:,:,i), stretchlim(reg_moving(:,:,i)), []);
        figure
        subplot(1,2,1)
        imshowpair(fixed_slice, moving_original_slice, 'Scaling', 'independent')
        title(sprintf('Slice %d before registration',i));
        subplot(1,2,2)
        imshowpair(fixed_slice, reg_moving_slice, 'Scaling', 'independent')
        title(sprintf('Slice %d after registration',i));
        
        figure
        imagesc(fixed(:,:,i))
        hold on
        b_tumour = bwboundaries(reg_roi(:,:,i));
        visboundaries(b_tumour,'Color','w','LineWidth', 1);
        %colormap(plasma)
        caxis([0,3500]);
        axis off
        daspect([1 1 1]);
        %colorbar(parula)
        title(sprintf('Slice %d with registered ROI',i));
    end
end



end
