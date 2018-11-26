function something = tnm034(im)

    % Image Capture
    % Read image file

    clc
    input_image = imread(im);

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
    figure, imshow(image_binary), hold on;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end

    %  Rotate

    most_frequent_angle = mode(T(P(:,2)));
    if(most_frequent_angle < 0)
        rotAngle = 90+most_frequent_angle;
    else
        rotAngle = most_frequent_angle - 90;
    end
    image_rotated = imrotate(image_binary, rotAngle);

    imshow(image_rotated);

    % Segmentation - Detect lines and save position
    % also - do some fine rotation adjustments 

    cropped_image = image_rotated(:,1:(size(image_rotated,2)/2));
    horizontal_summation = sum(cropped_image, 2);
    [pks, locs] = findpeaks(horizontal_summation);
    plot(pks , locs, '-o')

    % rotate by 0.01 degrees to find max
    highest_maxes = 0;
    rot_degree = 0;

    for degree = -0.1:0.01:0.1
        temp_sum = sum(imrotate(cropped_image, degree), 2);
        [temp_pks, temp_locs] = findpeaks(temp_sum);
        highest_val = sum(maxk(temp_pks, 5));
        if(highest_val > highest_maxes)
            highest_peaks = temp_pks;
            highest_locs = temp_locs;
            rot_degree = degree;
            highest_maxes = highest_val;
        end
    end

    new_rotated_image = imrotate(cropped_image, rot_degree);
    %figure, imshow(new_rotated_image);

    % save positions
    peak_threshold = floor(0.7* max(highest_peaks));
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

    staff = sortrows(staff);

    something = size(staff, 1);

    %% Segmentation - Thresholding



    %% Segmentation - Cleaning up (removing false objects)

    %% Segmentation - Labeling



    %% Classification
    % - Criteria:  length, area, circumference, form, color, texture
    % - Decision theory

    %% Symbolic description (output)

