function corners = get_corners(im_mask)
%Get the corners of a binary plate mask.

    % get plate corners via Harris Features
    corner_points = detectHarrisFeatures(im_mask);
    corner_points = corner_points.selectStrongest(10).Location;
    
    % eliminate any false points- in this case, eliminate points outside 
    % of the top/bottom 10% of columns and rows
    col_range = max(corner_points(:,1))-min(corner_points(:,1));
    row_range = max(corner_points(:,2))-min(corner_points(:,2));
    upper_col = max(corner_points(:,1))-col_range*.1;
    lower_col = min(corner_points(:,1))+col_range*.1;
    upper_row = max(corner_points(:,2))-row_range*.1;
    lower_row = min(corner_points(:,2))+row_range*.1;
    
    corner_points_outlier = (corner_points(:,2)<=lower_row | ...
        corner_points(:,2)>=upper_row) & ...
        (corner_points(:,1)<=lower_col | ...
        corner_points(:,1)>=upper_col);
    
    corner_points(~corner_points_outlier,:) = [];

    % classify each point
    [r,c] = find(im_mask);
    row_center = mean(r);
    col_center = mean(c);
    tl = []; tr = []; bl = []; br = [];
    for k = 1:length(corner_points)   
        % case 1: top left
        if corner_points(k,1) < row_center && corner_points(k,2) < col_center
            if isempty(tl)
                tl = corner_points(k,:);
            elseif corner_points(k,1) < tl(1) && corner_points(k,2) < tl(2)
                tl = corner_points(k,:);
            end
        % case 2: top right  
        elseif corner_points(k,1) > row_center && corner_points(k,2) < col_center
            if isempty(tr)
                tr = corner_points(k,:);
            elseif corner_points(k,1) > tr(1) && corner_points(k,2) < tr(2)
                tr = corner_points(k,:);
            end
        % case 3: bottom right
        elseif corner_points(k,1) > row_center && corner_points(k,2) > col_center
            if isempty(br)
                br = corner_points(k,:);
            elseif corner_points(k,1) > br(1) && corner_points(k,2) > br(2)
                br = corner_points(k,:);
            end
        % case 4: bottom left
        else
            if isempty(bl)
                bl = corner_points(k,:);
            elseif corner_points(k,1) < bl(1) && corner_points(k,2) > bl(2)
                bl = corner_points(k,:);
            end
        end
    end
    
    %return corner values
    corners = [tl; tr; br; bl];
    
end

