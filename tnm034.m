function outstr = tnm034(im)
 
    %input_image = imread('./Images_Training/im6s.jpg');
    input_image = imread(im);
    
    lookup_note_table = [ 'e4'; 'd4'; 'c4'; 'b3'; 'a3'; 'g3';
                          'f3'; 'e3'; 'd3'; 'c3'; 'b2'; 'a2';
                          'g2'; 'f2'; 'e2'; 'd2'; 'c2'; 'b1';
                          'a1'; 'g1'];
             
    output_chars = '';
    
    % Convert to grayscale and make binary
    image_binary = imcomplement(rgb2gray(input_image));
    image_binary = imbinarize(image_binary, 'adaptive');

    % Thicken to make thin lines thicker
    thickened_image = bwmorph(image_binary,'thicken', 1) ;

    % Use Hough to detect straight lines
    [H,T,~] = hough(thickened_image,'Theta',-90: 0.1 : 89.9);
    P  = houghpeaks(H, 1000); 

    %  Rough rotation - based on most frequent angles detected by Hough
    most_frequent_angle = mode(T(P(:,2)));
    if(most_frequent_angle < 0)
        rotAngle = 90+most_frequent_angle;
    else
        rotAngle = most_frequent_angle - 90;
    end
    image_rotated = imrotate(image_binary, rotAngle, 'bicubic');

    % Fine rotation - fine tune the rotation by rotating 0.1 degrees in
    % both direction until peaks are found, using horizontal projection on
    % half of the image (consider last  staff can be short)
    cropped_image = image_rotated(:,1:floor((size(image_rotated,2)/2)));
    highest_maxes = 0;
    rot_degree = 0;
    for degree = -0.1:0.01:0.1
        temp_sum = sum(imrotate(cropped_image, degree, 'bicubic'), 2);
        [temp_pks, temp_locs] = findpeaks(temp_sum);     
        sv = sort(temp_pks, 'descend');
        highest_val = sum(sv(1:5));
        if(highest_val > highest_maxes)
            highest_peaks = temp_pks;
            highest_locs = temp_locs;
            rot_degree = degree;
            highest_maxes = highest_val;
        end
    end

    % Save positions of peaks
    sorted_vals = sort(highest_peaks,'descend');
    top_ten = sorted_vals(1:10);
    
    peak_threshold = floor(0.7* median(top_ten));
    % thresholds the values so that the lower ones are zero
    % to find the correct staff lines and remove unwanted stuff
    filtered_peaks = highest_peaks.*(highest_peaks>peak_threshold);
    y_positions = highest_locs(filtered_peaks>0);

    % Some variables needed for sorting out staff lines
    no_rows = 1;
    staff = [];
    peak_size = size(y_positions,1);
    tolerable_distance = 5000;

    % Looping through all our maxima, to extract positions of the lines
    while(peak_size>=5)
        cluster = y_positions(1:5);
        previous_distance = cluster(5,1) - cluster(1,1);
        for pos = 1:1:peak_size-5
            % Here we check if the distance between the first found line and
            % the fifth found line in the cluster
            if(previous_distance > (y_positions(5+pos) - y_positions(1+pos)) && (y_positions(5+pos) - y_positions(1+pos)) < tolerable_distance)
                cluster = y_positions(1+pos:5+pos);
                previous_distance = (y_positions(5+pos) - y_positions(1+pos));
            end
        end
        staff(no_rows,:) = cluster(:);
        no_rows = no_rows + 1;
        % Deleting rows from the y_positions
        for ind = 1:1:5
            y_positions(y_positions == cluster(ind,1)) = [];
        end
        tolerable_distance = floor(previous_distance*1.2);
        peak_size = size(y_positions,1);
    end

    staffs = sortrows(staff);

    % rotate the original image for optimal precision
    new_rotated_image = imrotate(input_image, rotAngle + rot_degree, 'bicubic');

    % extract subimages, one for each staff, and scale them to have same
    % height, in order for all subimages to have the same relative staff line
    % positons
    margin_up = 0.88;
    margin_down = 0.75;
    img_height = size(new_rotated_image,1);
    binarize_threshold = 0.35;
    size_sub_image = 120;
    
    for staff_no = 1:size(staffs, 1)
        % Extract each staff
        distance = (staffs(staff_no, 5) - staffs(staff_no, 1) + 1);
        padding_top = 0;
        padding_bot = 0;
        
        start_pixel =staffs(staff_no, 1) - round(distance * margin_up);
        end_pixel = staffs(staff_no, 5) +  round(distance * margin_down);
        
        if (start_pixel < 1)
            padding_top = 1-start_pixel;
            start_pixel = 1;
        elseif (end_pixel > img_height)
            padding_bot = end_pixel - img_height;
            end_pixel = img_height;
        end
        
        % Pad image to avoid getting out of bounds on images where the top
        % and bottom of image is close to the edge of the image
        sub_image = new_rotated_image(start_pixel:end_pixel, :,:);
        sub_image = padarray(sub_image,[padding_top 0], 255, 'pre');
        sub_image = padarray(sub_image,[padding_bot 0], 255, 'post');
             
        % Resize image to a proper size
        factor = size_sub_image /size(sub_image,1);
        scaled = imresize(sub_image,factor,'bicubic');
        image_grayscale = imcomplement(rgb2gray(scaled));
        
        % find matching note heads using template matching. the template
        % comes from the training data.
        note_head_template = imcomplement(rgb2gray(imread('./note_head_9.png')));
        correlation = normxcorr2(note_head_template, image_grayscale);
        filtered_correlation = correlation > 0.5;
        
        % move results to center of note head, since the correlation
        % results are shifted
        filtered_correlation = circshift(filtered_correlation, [-round(size(note_head_template,1)/2), -round(size(note_head_template,2)/2)]);
        image_binary = imbinarize(image_grayscale, binarize_threshold);

        % Shift staff line coordinate system
        new_staff = staffs(staff_no, :)-start_pixel;
        new_staff = round(new_staff.*factor);      
        
        % Remove staff lines by opening with vertical kernel
        kernel_matrix = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
        opened_image = imopen(image_binary,kernel_matrix); 
        % Extract Area and Centroid properties from bwlabel
        areas = cell2mat(struct2cell(regionprops(filtered_correlation,'Area')))';
        centroids = regionprops(filtered_correlation,'centroid');
        sidemargin = size(note_head_template,1);  
        % Filter out outliers (non-notes) using area
        areas = areas.*(areas<25);
        area_threshold = 0.5 * max(areas);
        % Create and extended staff line vector that includes an element
        % for each staff line , ranging from e4 to g1
        staff_step = (new_staff(5)-new_staff(1))/8;
        extended_staff = [];
        extended_staff(1:6) = (new_staff(1)-(staff_step*6)):staff_step:new_staff(1)-staff_step;
        extended_staff(7) = new_staff(1);
        extended_staff(8) = new_staff(1) + ((new_staff(2)-new_staff(1))/2);
        extended_staff(9) = new_staff(2);
        extended_staff(10) = new_staff(2) + ((new_staff(3)-new_staff(2))/2);
        extended_staff(11) = new_staff(3);
        extended_staff(12) = new_staff(3) + ((new_staff(4)-new_staff(3))/2);
        extended_staff(13) = new_staff(4);
        extended_staff(14) = new_staff(4) + ((new_staff(5)-new_staff(4))/2);
        extended_staff(15) = new_staff(5);
        extended_staff(16:20) = new_staff(5)+staff_step:staff_step:new_staff(5)+(staff_step*5);
        for c = 1:size(centroids,1)
            % centroid position (middle of note head)
            position = round(cell2mat(struct2cell(centroids(c)))');    
            if(areas(c) < area_threshold) %don't make subimage if too small
                continue
            end
            
            % make a sub image for each note
            left = max(1,(position(1)-sidemargin));
            right = min((position(1)+sidemargin), size(opened_image,2));
            subimage = opened_image(:,left:right);
            % remove stem by opening
            horizontal_kernel = [ 0 0 0 0 0 ; 1 1 1 1 1; 0 0 0 0 0];
            subimage = imopen(subimage,horizontal_kernel);

            % make horizontal projection by summation
            horiz_proj_subimg = sum(subimage, 2);
            subarray = horiz_proj_subimg(position(2)-8:position(2)+8);
            % check if it is a note head
            centroid_width = horiz_proj_subimg(position(2));
            % don't classify if the detected centroid is not locatedin a note head 
            % ( it is too small or too big to be a note head)
            if(centroid_width < 9 || centroid_width > 17 || (min(subarray) > 0))
                continue; 
            end
         
            % check note head position to determine pitch
            [~, index] = min(abs(extended_staff-position(2)));   
                    
            % remove peaks below a certain value
            filter = horiz_proj_subimg > 6;
            horiz_proj_subimg = horiz_proj_subimg.*filter;
            
            % Convolve in order to "smoothen" spikey peaks
            conv_filter = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 ];
            conv_peaks = conv(horiz_proj_subimg, conv_filter);         
            [pks, ~] = findpeaks(conv_peaks);
            
             % if there is exactly one peak, return quarter, if two peaks 
             % return eigth
             if(length(pks) == 1)
                 output_chars = strcat(output_chars, upper(lookup_note_table(index,:)));
             elseif(length(pks) == 2)
                 output_chars = strcat(output_chars, lookup_note_table(index,:));         
             end           
        end
        output_chars = strcat(output_chars, 'n');
    end

    outstr = convertCharsToStrings(output_chars)
