%function something = tnm034(im)
    
    % Image Capture
    % Read image file
    clear
    clc
    input_image = imread('./Images_Training/Le_1_Example.jpg');
    %input_image = imread(im);

    % Preprocessing - Make binary

    % convert to grayscale and invert
    image_binary = imcomplement(rgb2gray(input_image));

    % binarizes image I with a global threshold computed using Otsu's method
    image_binary = imbinarize(image_binary, 'adaptive');
    % testa andra trösklar, ljusgråa linjer kommer med nu

    % thicken to make thin lines really POP
    thickened_image = bwmorph(image_binary,'thicken', 1) ;

    % Preprocessing
    % - Geometric restoration,Photometric  restoration,Noise removal?

    % Preprocessing - Detection of lines and edges, to rotate
    % - Staff lines must be located and rotated. May be removed.

    [H,T,R] = hough(thickened_image,'Theta',-90: 0.1 : 89.9);
    P  = houghpeaks(H, 1000); 
    %imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    %xlabel('\theta'), ylabel('\rho');
    %axis on, axis normal, hold on;
    % show white boxes where peaks are
    %plot(T(P(:,2)),R(P(:,1)),'s','color','white');

    %
    % min length of line is based on image width
    min_length_staff_lines = floor(0.2*size(input_image,2));

    lines = houghlines(thickened_image,T,R,P,'FillGap',10,'MinLength',min_length_staff_lines);
    %figure, imshow(image_binary), hold on;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end

    %  Rotate

    most_frequent_angle = mode(T(P(:,2)));
    if(most_frequent_angle < 0)
        rotAngle = 90+most_frequent_angle;
    else
        rotAngle = most_frequent_angle - 90;
    end
    image_rotated = imrotate(image_binary, rotAngle, 'bicubic');

    % Segmentation - Detect lines and save position
    % also - do some fine rotation adjustments 
    
    cropped_image = image_rotated(:,1:floor((size(image_rotated,2)/2)));
    horizontal_summation = sum(cropped_image, 2);
    [pks, locs] = findpeaks(horizontal_summation);
    %plot(pks , locs, '-o')

    % rotate by 0.01 degrees to find max
    highest_maxes = 0;
    rot_degree = 0;

    for degree = -0.1:0.01:0.1
        temp_sum = sum(imrotate(cropped_image, degree, 'bicubic'), 2);
        [temp_pks, temp_locs] = findpeaks(temp_sum);
        
        % For premium users only!
        %highest_val = sum(maxk(temp_pks, 5)); F
        sv = sort(temp_pks, 'descend');
        highest_val = sum(sv(1:5));
        if(highest_val > highest_maxes)
            highest_peaks = temp_pks;
            highest_locs = temp_locs;
            rot_degree = degree;
            highest_maxes = highest_val;
        end
    end

    %new_rotated_cropped_image = imrotate(cropped_image, rot_degree);
    %new_rotated_image = imrotate(image_rotated, rot_degree);

    % save positions
    % Premium user? Use maxk!
    sorted_vals = sort(highest_peaks,'descend');
    top_ten = sorted_vals(1:10);
    
    peak_threshold = floor(0.7* median(top_ten));
    % thresholds the values so that the lower ones are zero
    filtered_peaks = highest_peaks.*(highest_peaks>peak_threshold);
    y_positions = highest_locs(filtered_peaks>0);

    % Some variables needed for sorting out staff lines
    no_rows = 1;
    staff = [];
    peak_size = size(y_positions,1);
    tolerable_distance = 5000;

    % Looping through all our maxima lines, to extract positions of the staff
    % lines
    while(peak_size>=5)
        cluster = y_positions(1:5);
        previous_distance = cluster(5,1) - cluster(1,1);
        for pos = 1:1:peak_size-5

            % Here we check if the distance between the first found line and
            % the last found line in the cluster
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

   % snurra original-bilden för bästa resultat
    new_rotated_image = imrotate(input_image, rotAngle + rot_degree, 'bicubic');

    % ta ut sub-bilder. en för varje staff. gör till samma storlek
    % så vi hela tiden får samma referenskoordinater, så kärnan kan ha samma storlek
    margin_up = 0.88;
    margin_down = 0.75;

    img_height = size(new_rotated_image,1);

    binarize_threshold = 0.35;
    size_sub_image = 120;
    
    counter = 0;
    
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
      
        sub_image = new_rotated_image(start_pixel:end_pixel, :,:);
        sub_image = padarray(sub_image,[padding_top 0], 255, 'pre');
        sub_image = padarray(sub_image,[padding_bot 0], 255, 'post');
             
        % Resize image to a proper size
        factor = size_sub_image /size(sub_image,1);
        scaled = imresize(sub_image,factor,'bicubic');
        %figure, imshow(scaled);     
        image_grayscale = imcomplement(rgb2gray(scaled));
        
        % find matching note heads
        %note_head_template = imcomplement(rgb2gray(imread('./note_head.png')));
        note_head_template = imcomplement(rgb2gray(imread('./note_head_9.png')));
        correlation = normxcorr2(note_head_template, image_grayscale);
        filtered_correlation = correlation > 0.5;
        filtered_correlation = circshift(filtered_correlation, [-round(size(note_head_template,1)/2), -round(size(note_head_template,2)/2)]);

  
        image_binary = imbinarize(image_grayscale, binarize_threshold);

        %figure, subplot(1,1,1), imshow(filtered_correlation) ,subplot(1,2,1), imshow(image_binary)
        
        
        
        new_staff = staffs(staff_no, :)-start_pixel;
        new_staff = round(new_staff.*factor);
        
        
        kernel_matrix = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
        opened_image = imopen(image_binary,kernel_matrix);
        %cleaned_image = bwmorph(opened_image, 'clean');
        %kernel_matrix = [0 0 0; 1 1 1; 0 0 0];
        %closed_image = imclose(opened_image,kernel_matrix);
        
        %horizontal_kernel = [ 0 0 0 0; 1 1 1 1; 0 0 0 0; 0 0 0 0];
        %opened_image = imopen(opened_image,horizontal_kernel);
        
        % IMAGE WITH ONLY NOTE LINES
        vertical_kernel = [0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0;
            0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0];
        %note_line_img = imopen(image_binary, vertical_kernel);
        %figure, imshow(note_line_img);
        
        % IMAGE WITH ONLY NOTE HEADS
        skewed_circle_kernel = [ 
            0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 1 1 1 1 1 0 0 0;
            0 0 0 1 1 1 1 1 1 1 0 0;
            0 0 1 1 1 1 1 1 1 1 1 0;
            0 1 1 1 1 1 1 1 1 1 1 1;
            1 1 1 1 1 1 1 1 1 1 1 0;
            0 1 1 1 1 1 1 1 1 1 0 0;
            0 0 1 1 1 1 1 1 1 0 0 0;
            0 0 0 1 1 1 1 1 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0;];
        %note_head_img = imopen(image_binary, skewed_circle_kernel);
        
        
        labeled_image = bwlabel(filtered_correlation);
        areas = cell2mat(struct2cell(regionprops(labeled_image,'Area')))';
        centroids = regionprops(labeled_image,'centroid');
        sidemargin = size(note_head_template,1);  
%          figure; imshowpair(image_grayscale,filtered_correlation);
%          hold on;
        outliers = areas<25;
        areas = areas.*outliers;
        area_threshold = 0.5 * max(areas);
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
            position = round(cell2mat(struct2cell(centroids(c)))');
            
            if(areas(c) < area_threshold) %don't make subimage if too small
                continue
            end
%             plot(position(1),position(2), 'or');
            
            left = max(1,(position(1)-sidemargin));
            right = min((position(1)+sidemargin), size(opened_image,2));
            subimage = opened_image(:,left:right);
            horizontal_kernel = [ 0 0 0 0 0 ; 1 1 1 1 1; 0 0 0 0 0];
            subimage = imopen(subimage,horizontal_kernel);
            %subimage = bwmorph(subimage, 'majority');

            vert_proj_subimg = sum(subimage, 2);
            [pks, locs] = findpeaks(vert_proj_subimg);
            %figure, plot(vert_proj_subimg);
            subarray = vert_proj_subimg(position(2)-8:position(2)+8);
            % check if it is a note head
            centroid_width = vert_proj_subimg(position(2));
            if(centroid_width < 9 || centroid_width > 17 || (min(subarray) > 0))

                continue;
            end
            
            % check note head position to determine pitch
            counter = counter + 1

            position(2)
            
            for i = 1:20
                pos = round(extended_staff(i));
                subimage(pos,:) = 1;
            end
            
            
            
            % remove peaks below a certain value
            filter = vert_proj_subimg > 6;
            vert_proj_subimg = vert_proj_subimg.*filter;
            
              figure('Name',num2str(counter)), imshow(subimage);
%              figure, plot(vert_proj_subimg);
            
                
        end
%          hold off; 

        
    end
    
    %figure, imshow(labeled_image);
    %figure, bar(vertical_summation);
    
    something = 'hej';
    
    %imshow(new_rotated_image);
    
    % Segmentation - Thresholding
    
    %% Segmentation - Cleaning up (removing false objects)

    %% Segmentation - Labeling



    %% Classification
    % - Criteria:  length, area, circumference, form, color, texture
    % - Decision theory

    %% Symbolic description (output)

