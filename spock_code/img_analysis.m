% High level strategy:
%   1. Cycle through entire library to collect population level stats
%   2. Loop through all images:
%       2.1 For each image, normalize to set brightness, lum, hue, etc.
%       2.2 Find all circles in each image, store stats about each
%       2.3 Store circles in standard format for later comparison
%   3. After all images processed, collect population level circle stats
%   4. Filter out any outliers from circle population
%   5. Compare cross-image for replicant similarity
%
% Notes:
%   - Set im_debug=1 to see image debugging for troubleshooting
%   - Script assumes data is all in {Current Directory}/Data
%   - All files are saved in {Current Directory}/Matfiles

clear all; close all;
im_debug = 1;
warning('off','images:imfindcircles:warnForLargeRadiusRange');
[~,~,gene_list] = xlsread('Collins_lab_Keio_map.xlsx');
welllocation = reshape([gene_list(2:385,3)].',[24,16]).';
fileList = dir('Data');

% Loop through all images to process
for k = 1:length(fileList)

    % check if current item is a directory
    if fileList(k).isdir
        continue
    elseif strcmp(fileList(k).name,'.DS_Store')
        continue
    end

    % check if it's a PC or Mac to build file structure
    if ispc
        im = imread(['Data\',fileList(k).name]);
    else
        im = imread(['Data/',fileList(k).name]);
    end

    % create mask of the plate
    % steps: create grayscale img, binarize, then close the circles inside
    % the binarized mask, and finally convert to uint16 for multiplication
    if k == 10 || k == 12
        im(:,:,4)=[];
    end
    im_gr = rgb2gray(im);
    im_mask = im2bw(im_gr, graythresh(im_gr));
    im_mask = imclose(im_mask,strel('disk',15));
    im_mask = uint16(im_mask);

    if im_debug
        close all;
        figure, imshow(imfuse(im_gr,im_mask));
        title('BW image with plate mask')
        hold on
    end

    % create contrasted image for easier processing
    % this is a contrast-adjusted image multiplied by the mask we created,
    % which gives a binary-zero outer area even in when complemented
    im_co = imcomplement(imadjust(im_gr.*im_mask)).*im_mask;

    % find circles in image; need multiple rounds to find both
    % low-intensity and high-intensity, try to keep rMax < 4*rMin
    [c1, r1] = imfindcircles(adapthisteq(im_co),[6 23]);
    [c2, r2] = imfindcircles(adapthisteq(im_gr,'clipLimit',0.03),[6 23]);

    % eliminate any repeated circles
    is_dup = false(size(c2,1),1);
    for j = 1:size(c2,1)
        is_dup(j) = any(in_circle(c1,r1,c2(j,1),c2(j,2)));
    end
    c2(is_dup,:) = [];
    r2(is_dup) = [];

    % combine centers and radii
    centers = [c1;c2];
    radii = [r1;r2];

    %create circle mask
    mask2 = zeros(size(im_mask));
    for j = 1:length(centers)
        xc = centers(j,1); yc = centers(j,2);
        xi = size(im,2); yi = size(im,1);
        rad = radii(j);
        [xx,yy] = meshgrid(1:xi,1:yi);
        mask = (xx - xc).^2 + (yy - yc).^2 <= rad.^2;
        mask2 = mask+mask2;
    end

    % image background and comp is found using colony mask (mask2)
    mean_bg = mean(im_gr(logical(im_mask.*uint16(~mask2))));
    experiment.background = mean_bg;
    mean_bg2 = mean(imcomplement(im_gr(logical(im_mask.*uint16(~mask2)))));
    experiment.background2 = mean_bg2;
    
    % second round of processing for circle-removed images
    im_gr2 = (im_gr-mean_bg).*uint16(~mask2);
    im_gr2 = imadjust(im_gr2,stretchlim(im_gr2,[0.0005 0.9995]),[]);
    im_co2 = imadjust((imcomplement(im_gr)-mean_bg2)).*im_mask.*uint16(~mask2);

    % find circles in image; need multiple rounds to find both
    % low-intensity and high-intensity, try to keep rMax < 4*rMin
    [c1, r1] = imfindcircles(im_co2,[6 23]);
    [c2, r2] = imfindcircles(im_gr2,[6 23]);
    c_new = [c1;c2];
    r_new = [r1;r2];

    is_dup = false(size(c_new,1),1);
    for j = 1:size(c_new,1)
        is_dup(j) = any(in_circle(centers,radii,c_new(j,1),c_new(j,2)));
    end
    c_new(is_dup,:) = [];
    r_new(is_dup) = [];

    centers = [centers; c_new];
    radii = [radii; r_new];

    % show image if flagged
    if im_debug
        figure, imshow(imfuse(im_co,im_mask));
        hold on
        viscircles(centers,radii,'EdgeColor','b');
        title('Contrast Image with Found Circles')
    end

    % set plate size
    num_rows = 16;
    num_cols = 24;
    num_wells = num_rows * num_cols;

	% get plate outline; click on bottom left and top left,
	% then bottom right to top right to outline shape
	if ~im_debug
		f = figure
		imshow(imfuse(im_co,im_mask));
	end
	
    pts_rd1 = getline;
    tl = pts_rd1(2,:);
    bl = pts_rd1(1,:);

    pts_rd2 = getline;
    tr = pts_rd2(2,:);
    br = pts_rd2(1,:);

    if im_debug
        plotPts = [tl;tr;br;bl;tl];
        plot(plotPts(:,1),plotPts(:,2),'g*-','MarkerSize',10,'LineWidth',2);
    else
		close(f)
	end

    % create grid structure
    % deltaX and deltaY are the slope needed to follow the plate
    col_deltaX = (tr(1)-tl(1))/num_cols;
    col_deltaY = (tr(2)-tl(2))/num_cols;
    row_deltaX = (bl(1)-tl(1))/num_rows;
    row_deltaY = (bl(2)-tl(2))/num_rows;
    grid = struct;

    % initialize first point
    grid(1,1).tl = tl;
    grid(1,1).tr = tl + [col_deltaX, col_deltaY];
    grid(1,1).br = grid(1,1).tr + [row_deltaX, row_deltaY];
    grid(1,1).bl = tl + [row_deltaX, row_deltaY];
    centers_tmp = centers;
    radii_tmp = radii;
    r=1; c=1;
    squarePts = [grid(r,c).tl;grid(r,c).tr;grid(r,c).br;grid(r,c).bl;grid(r,c).tl];
    if im_debug
        plot(squarePts(:,1),squarePts(:,2),'r');
    end
    in_sq = inpolygon(centers_tmp(:,1),centers_tmp(:,2),squarePts(:,1),squarePts(:,2));
    grid(r,c).centers = centers_tmp(in_sq,:);
    grid(r,c).radii = radii_tmp(in_sq);
    centers_tmp(in_sq,:) = [];
    radii_tmp(in_sq) = [];

    % build grid structure
    for r = 1:num_rows
        for c = 1:num_cols
            if c==1 && r == 1
                continue
            elseif c == 1
                grid(r,c).tl = grid(r-1,c).bl;
            else
                grid(r,c).tl = grid(r,c-1).tr;
            end
            grid(r,c).tr = grid(r,c).tl + [col_deltaX, col_deltaY];
            grid(r,c).br = grid(r,c).tr + [row_deltaX, row_deltaY];
            grid(r,c).bl = grid(r,c).tl + [row_deltaX, row_deltaY];
            squarePts = [grid(r,c).tl;grid(r,c).tr;grid(r,c).br;grid(r,c).bl;grid(r,c).tl];

            if im_debug
                plot(squarePts(:,1),squarePts(:,2),'r');
            end

            % find any circles inside the current grid square
            % only the center of the circle must fall inside the grid!
            in_sq = inpolygon(centers_tmp(:,1),centers_tmp(:,2),squarePts(:,1),squarePts(:,2));
            grid(r,c).centers = centers_tmp(in_sq,:);
            grid(r,c).radii = radii_tmp(in_sq);
            centers_tmp(in_sq,:) = [];
            radii_tmp(in_sq) = [];
        end
    end

    experiment.grid = grid;
    experiment.image = im;
    experiment.image_mask = im_mask;
    experiment.image_name = fileList(k).name;

    % save current grid as the drug and plate name
    file_split = strsplit(fileList(k).name,' ');
    file_split2 = strsplit(file_split{4},'.');
    drugname = file_split{1};
    drugnum = file_split{2};
    platenum = file_split2{1};
    saveName = [drugname,'_',drugnum,'_plate_',platenum];
    eval([saveName,'= experiment;']);
    if ispc
        eval(['save(''',pwd,'\Matfiles\',saveName,''',''',saveName,''')']);
    else
        eval(['save(''',pwd,'/Matfiles/',saveName,''',''',saveName,''')']);
    end

    if im_debug
        input('Press enter to continue...');
    end
end


% Cycle through entire all mat files to collect stats on colonies
fileList = dir('Matfiles');
plate_intensity_sub = zeros(1,length(fileList));
background_intensity = zeros(1,length(fileList));
for k = 1:length(fileList)
    
    % check if current item is a directory
    if fileList(k).isdir
        continue
    elseif strcmp(fileList(k).name,'.DS_Store')
        continue
    end
    
    % check if it's a PC or Mac to build file structure
    if ispc
        str = load([fileList(k).folder,'\',fileList(k).name]);
    else
        str = load([fileList(k).folder,'/',fileList(k).name]);
    end
    
    %dynamically load struct from file
    structName = fieldnames(str);
    structName = structName{1};
    experiment = str.(structName);
    experiment.grid = (experiment.grid); %flipud(experiment.grid)
    I = rgb2gray(experiment.image);
    experiment.I = I;
    grid = experiment.grid;
    
    all_intensities = [];
    for r = 1:size(grid,1)
        for c = 1:size(grid,2)
            
            % take only colony closest to center of grid
            if length(grid(r,c).radii)>1
                dist=[];
                x1 = grid(r,c).tl(1); x2 = grid(r,c).br(1);
                y1 = grid(r,c).tl(2); y2 = grid(r,c).br(2);
                M = [(x1+x2)./2,(y1+y2)./2];
                for z = 1:length(grid(r,c).radii)
                    x2test = grid(r,c).centers(z,1);
                    y2test = grid(r,c).centers(z,2);
                    dist(z) = sqrt((x2test-M(1)).^2+(y2test-M(2))^2);
                end
                [xx ind1] = (min(dist));
                [xx ind2] = (max(grid(r,c).radii));
                if ind1 == ind2
                    ind = ind2;
                else
                    ind = ind1;
                end
                
                grid(r,c).radii = grid(r,c).radii(ind);
                grid(r,c).centers = grid(r,c).centers(ind,:);
                
            elseif isempty(grid(r,c).radii)
                grid(r,c).centers = [0 0];
                grid(r,c).radii = 0;
            end
            
            % build statistics and write to struct
            xc = grid(r,c).centers(1); yc = grid(r,c).centers(2);
            xi = size(experiment.image,1); yi = size(experiment.image,2);
            rad = grid(r,c).radii;
            [xx,yy] = meshgrid(1:yi,1:xi);
            mask = (xx - xc).^2 + (yy - yc).^2 <= rad.^2;
            intensities = I(mask);
            avg_intensity = mean(intensities);
            if isnan(avg_intensity); avg_intensity = 0; end
            stdavg_intensity = std(double(intensities));
            experiment.grid(r,c).mean_colony_intensity = avg_intensity;
            experiment.grid(r,c).std_colony_intensity = stdavg_intensity;
            experiment.I = I;
            all_intensities(end+1) = mean(intensities);
        end
    end
    
    all_intensities(isnan(all_intensities)) = 0;
    experiment.plate_intensity_mean = mean(all_intensities);
    eval(['save(''',pwd,'/Matfiles/',fileList(k).name,''',''','experiment',''')']);
end

% Cycle through all mat files to filter out invalid colonies
fileList = dir('Matfiles');
colavg = [];
for k = 1:length(fileList)
    
    % check if current item is a directory
    if fileList(k).isdir
        continue
    elseif strcmp(fileList(k).name,'.DS_Store')
        continue
    end
    
    % check if it's a PC or Mac to build file structure
    if ispc
        str = load([fileList(k).folder,'\',fileList(k).name]);
    else
        str = load([fileList(k).folder,'/',fileList(k).name]);
    end
    
    %dynamically load struct from file
    structName = fieldnames(str);
    structName = structName{1};
    experiment = str.(structName);
    grid = experiment.grid;
    
    plate_intensity = reshape([experiment.grid(:,:).mean_colony_intensity],[16,24]);
    plate_intensity_ind = plate_intensity~=0;
    
    % build row/col means
    rowmeans = zeros(size(grid,1),size(grid,2));
    colmeans = zeros(size(grid,1),size(grid,2));
    for r = 1:size(grid,1)
        for c = 1:size(grid,2)
            indR = plate_intensity_ind(r,:);
            indC = plate_intensity_ind(:,c);
            pR = plate_intensity(r,:);
            pC = plate_intensity(:,c);
            pR = pR(indR); pC = pC(indC);
            rowmeans(r,c) = mean(pR);
            colmeans(r,c) = mean(pC);
        end
    end
    
    % normalize intensity by row/col
    plotintensities = zeros(size(grid,1),size(grid,2));
    for r = 1:size(grid,1)
        for c = 1:size(grid,2)
            if plate_intensity(r,c)==0
                continue
            end
            fluor_per_cell = experiment.grid(r,c).mean_colony_intensity;%./mean(mean(plate_intensity(plate_intensity>0)));
            m1(r,c) = fluor_per_cell;
            fluor_per_cell = fluor_per_cell./rowmeans(r,c);
            m2(r,c) = fluor_per_cell;
            fluor_per_cell = fluor_per_cell./colmeans(r,c);
            m3(r,c) = fluor_per_cell;
            
            if isnan(fluor_per_cell); fluor_per_cell = 0; end
            experiment.grid(r,c).mean_colony_intensity = fluor_per_cell;
            plotintensities(r,c) = fluor_per_cell;
        end
    end
    
    % update struct fields
    totalaverage = mean(mean(plotintensities(plotintensities>0)));
    for r = 1:size(grid,1)
        for c = 1:size(grid,2)
            fluor_per_cell = experiment.grid(r,c).mean_colony_intensity./totalaverage;
            if isnan(fluor_per_cell); fluor_per_cell = 0; end
            experiment.grid(r,c).mean_colony_intensity = fluor_per_cell;
            plotintensities(r,c) = fluor_per_cell;
            colavg(end+1) = fluor_per_cell;
        end
    end
	
    if im_debug
        figure
        subplot(2,2,1);     imagesc(m1)
        subplot(2,2,2);     imagesc(m2)
        subplot(2,2,3);     imagesc(m3)
        subplot(2,2,4);     imagesc(plotintensities)
    end
    
    % write gene info to struct
    num_cells = 16*24;
    start_idx = (k-3)*num_cells+2;
    end_idx = (k-2)*num_cells+1;
    curr_genePlate = reshape([gene_list(start_idx:end_idx,4)].',[24,16]).';
    curr_geneName = reshape([gene_list(start_idx:end_idx,2)].',[24,16]).';
    for r = 1:size(grid,1)
        for c = 1:size(grid,2)
            experiment.grid(r,c).geneLoc = welllocation(r,c);
            experiment.grid(r,c).genePlate = curr_genePlate(r,c);
            experiment.grid(r,c).geneName = curr_geneName(r,c);
        end
    end
    eval(['save(''',pwd,'/Matfiles/',fileList(k).name,''',''','experiment',''')']);
    
end

save colavg1 colavg
figure; plot(colavg,'o')
figure; hist(colavg,50), set(gca,'fontsize',24),title('colony pixel distribution')
