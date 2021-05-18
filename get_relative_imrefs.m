function [R1, R2] = get_relative_imrefs(info1, info2, reg_type)
%% Get relative spatial references for the two images to be registered

% Ambre's edits from the original script: 
% modified to allow for spatial references of 4D timeseries data to be extracted
% reg_type to perform 2D to 2D reg on 3D data if slice extraction is required


info = [info1, info2];
dimensions = [length(info1.ImageSize), length(info2.ImageSize)];

% reduce dimensionality 
for i=1:2
    if dimensions(i) == 4
        dimensions(i) = 3;   % in case the image is a 4D timeseries
    end
    if reg_type == 2
        dimensions(i) = 2;
    end
end


% Initialise Q that will store the coordinates for four corners of the
% image
Q = zeros(3, 4*length(info));

% Assuming images have the same dimensionality, either info1 or info2 can
% be used in next line
if dimensions(1) == 3
    Q = [Q; zeros(1, 4*length(info))];
end


%% For 3D registration there must be a minimum of 16 slices
for i = 1:length(info)
    padding = 0;
    if dimensions(1) == 3
        z_pixels = info(i).ImageSize(3);
        if info(i).ImageSize(3) < 16
            padding = round((16 - info(i).ImageSize(3))/2);
            info(i).ImageSize(3) = info(i).ImageSize(3) + 2*padding;
        end
    end
    
    
   %% Calculate world coordinates by applying the transformations stored in the header
    x_pixels = info(i).ImageSize(1);
    y_pixels = info(i).ImageSize(2);
       
    P = [-0.5, -0.5, -0.5-padding 1;                       %top left corner of bottom layer
        x_pixels-0.5, -0.5, -0.5-padding 1;               %top right corner of bottom layer
        -0.5, y_pixels-0.5, -0.5-padding 1];               %bottom left corner of bottom layer
        
    if dimensions(1) == 3
        P = [P; -0.5, -0.5, z_pixels-0.5+padding 1];              %top left corner of top layer
    end
    
    %Location of corners in world coordinates
    Q(:,1+(i-1)*4:i*4) = P * info(i).Transform.T;
end

%Remove the columns that were only needed to apply affine transformation
Q = Q(:, [1:3, 5:7]);

%unit vectors in the original x y and z directions after the transformation
n_x = (Q(2, 1:3)-Q(1, 1:3))/norm(Q(2, 1:3)-Q(1, 1:3));
n_y = (Q(3, 1:3)-Q(1, 1:3))/norm(Q(3, 1:3)-Q(1, 1:3));

%Determine the Image Extent In World Coordinates
x1 = norm(Q(2, 1:3)-Q(1, 1:3));
y1 = norm(Q(3, 1:3)-Q(1, 1:3));
x2 = norm(Q(2, 4:6)-Q(1, 4:6));
y2 = norm(Q(3, 4:6)-Q(1, 4:6));

%Determine the relative distance between the upper left corners of the
%two resulting images
dx = dot(n_x, Q(1, 4:6)-Q(1, 1:3));
dy = dot(n_y, Q(1, 4:6)-Q(1, 1:3));

if dimensions(2) == 3
    n_z = (Q(4, 1:3)-Q(1, 1:3))/norm(Q(4, 1:3)-Q(1, 1:3));
    z1 = norm(Q(4, 1:3)-Q(1, 1:3));
    z2 = norm(Q(4, 4:6)-Q(1, 4:6));
    dz = dot(n_z, Q(1, 4:6)-Q(1, 1:3));
    R1 = imref3d(info(1).ImageSize(1:3), [0 x1], [0 y1], [0 z1]);
    R2 = imref3d(info(2).ImageSize(1:3), [dx x2+dx], [dy y2+dy], [dz z2+dz]);
else
    R1 = imref2d(info(1).ImageSize(1:3), [0 x1], [0 y1]);
    R2 = imref2d(info(2).ImageSize(1:3), [dx x2+dx], [dy y2+dy]);

end



end