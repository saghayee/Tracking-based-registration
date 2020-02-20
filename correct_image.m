function IMG_reg = correct_image(IMG,a,pad_d)


%         im_red = imread([savefolder '\redChan\red' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
%         im_green = imread([savefolder '\greenChan\green' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %end
        % Prior to translation, pad both dimensions with random background,
    % with a pad width equal to the maximum displacement; after using
    % imtranslate, cut this padded region out
    
    % Currently, padding just repeats the the border onto the added rows
    % and columns; the size of the padding is determined by the maximum
    % size in translational drift + the maximum drift that can result from
    % rotation
    IMG = padarray(IMG,pad_d*ones(1,2),'replicate');
    
    %im_red = padarray(im_red,pad_d*ones(1,2),'replicate');
    %im_green = padarray(im_green,pad_d*ones(1,2),'replicate');
    IMG_trans = imtranslate(IMG,[a(1) a(2)],'nearest');
    %im_red_trans = imtranslate(im_red,[a_all(i,1) a_all(i,2)]);
    %im_green_trans = imtranslate(im_green,[a_all(i,1) a_all(i,2)]);
    % In order to properly rotate image, the translated image must be
    % padded so that the cells' centroid corresponds to the image center
    if a(3) ~= 0
        IMG_reg_sl = rotateAround(IMG_trans,a(5)+pad_d,a(4)+pad_d,180/pi*a(3),'nearest');
        
        %im_red_reg = rotateAround(im_red_trans,a_all(i,8)+pad_d,a_all(i,7)+pad_d,180/pi*a_all(i,3)*sign(a_all(i,6)),'bilinear');
        %im_green_reg = rotateAround(im_green_trans,a_all(i,8)+pad_d,a_all(i,7)+pad_d,180/pi*a_all(i,3)*sign(a_all(i,6)),'bilinear');
    else
        IMG_reg_sl = IMG_trans;
        
        %im_red_reg = im_red_trans;
        %im_green_reg = im_green_trans;
    end
    IMG_reg_sl(1:pad_d,:) = [];
    IMG_reg_sl(:,1:pad_d) = [];
    IMG_reg_sl(end-pad_d+1:end,:) = [];
    IMG_reg_sl(:,end-pad_d+1:end) = [];
    IMG_reg = IMG_reg_sl;