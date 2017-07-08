clc; clf; clear all; close all;

% Define global variables used for constant settings
global showLabels showText testRun showFigures dataDir templateChannel targetChannel thresholdChannel maxThreshold minThreshold

% Channels used for analysis
templateChannel = 1;
targetChannel = 2;
thresholdChannel = 0;
% Threshold values (1 = max, 0 = min)
maxThreshold = 1.0;
minThreshold = 0.05;
% Show (1) or hide (0) figures during analysis
showFigures = 1; % Show analysis images
showText = 1; % Show cell label info
showLabels = 1; % Show cell labels
testRun = 1; % Do a testrun of the first image in the first folder
dataDir = 'data'; % Folder containing data files

% Adjust for matlabs index 1
templateChannel = templateChannel+1;
targetChannel = targetChannel+1;
thresholdChannel = thresholdChannel+1;
% Adjust threshold based on 16-bit depth
maxThreshold = maxThreshold*2^16;
minThreshold = minThreshold*2^16;

set(0, 'DefaulttextInterpreter', 'none') % Disable LaTeX syntax in figure texts

% Get all subdirs in the data dir
files = dir(dataDir);
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

% Loop through the subdirs and analyze the files
for i = 1:length(directoryNames)
    analyzefiles(directoryNames{i});
    if testRun
        break
    end
end

% Analyzes all .czi files in a directory
function analyzefiles(directory)
    global testRun dataDir

    files = dir(strcat(dataDir, '/', directory, '/*.czi'));
    cells = [];
    for i = 1:length(files)
        file = files(i).name;
        folder = files(i).folder;
        data = bfopen(strcat(folder, '/', file)); % Opens the .czi file
        
        % Gets the um/pixel scale from czi metadata
        metadata = data{1, 4};
        scale = metadata.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
        
        series = data{1, 1}; % Gets the data set from the .czi file
        % Gets the average brightness of each cell, measured in each of the
        % channels (template, target, threshold) as a N by 3 matrix, where
        % N is the number of cells found in the template channel. The cells
        % brightness is merged with the other cells found from the .czi
        % files in the same folder
        cells = [cells; cellcountfluorescence(series, file, scale)];
        if testRun
            break
        end
    end
    
    if(~isempty(cells)) % Only do ratio calculations if cells were found
        % Gets the FRET as a ratio between individual cell brightness
        ratio = cells(:,2)./(cells(:,1)+cells(:,2));
        average = mean(ratio);
        deviation = std(ratio);

        % Show FRET calculation, and save figure and data
        figure
        histogram(ratio, 0:0.05:1);
        title(directory);
        string = strcat(num2str(average), '±', num2str(deviation));
        text(0.5,0.95,{'Average ratio', string},'Units','normalized','HorizontalAlignment', 'center');
        saveas(gcf,strcat('results/', directory, '_histogram.png'));
        csvwrite(strcat('results/', directory, '.csv'), cells);
    end
end

% Find the cells brightness using a fluorescence channel as template for
% identifying cells
function [cells] = cellcountfluorescence(series, filename, scale)
    global templateChannel targetChannel thresholdChannel

    targetImage = series{targetChannel, 1};
    templateImage = series{templateChannel, 1};
    thresholdImage = series{thresholdChannel, 1};

    im_bw = im2bw(templateImage, 0.1); % Converts to binary image ignoring brightness below threshold
    im_open = imopen(im_bw, strel('disk',1)); % Removes single pixel noise
    im_open = bwareaopen(im_open, 1); % Removes cells below pixel size
    
    % Uses watershed segmentation to separate cells that are clumped
    % together
    im_max = imextendedmax(templateImage, 12000);
    im_dist = -bwdist(~im_open); 
    im_dist = imimposemin(im_dist, im_max);
    im_watershed = watershed(im_dist);
    im_watershedadjusted = im_open;
    im_watershedadjusted(im_watershed == 0) = 0;

    % Labels all the cells found in the watershed segmented image
    [im_labels, label_count] = bwlabel(im_watershedadjusted, 4);

    % Gets the average brightness of each of the cells in the channels
    [cells] = getbrightness(im_labels, label_count, templateImage, targetImage, thresholdImage, filename, scale);

end

% Gets the average brightness of each of the cells in the channels
function [cells] = getbrightness(im_labels, count, templateImage, targetImage, thresholdImage, filename, scale)
    global showLabels showFigures minThreshold maxThreshold

    cells = [];
    newregions = zeros(512, 512); % Placeholder for the cells that passes threshold selection
    oldregions = zeros(512, 512); % Placeholder for the cells that doesn't pass threshold selection
    
    % Loop through all the labels found in the segmented cell image
    for i = 1:count
        % Gets the pixels in the label
        [r, c] = find(im_labels==i);
        rc = [r c];
        pixels = size(rc, 1);
        sumTemplate = 0;
        sumTarget = 0;
        sumThreshold = 0;
        fret = 0;

        % Loop through all the pixels in the label, and add their intensity
        for j = 1:pixels
            valueTemplate = double(templateImage(rc(j,1), rc(j,2)));
            valueTarget = double(targetImage(rc(j,1), rc(j,2)));
            valueThreshold = double(thresholdImage(rc(j,1), rc(j,2)));
            
            sumTemplate = sumTemplate + valueTemplate;
            sumTarget = sumTarget + valueTarget;
            sumThreshold = sumThreshold + valueThreshold;
        end
        
        % Divide by the pixel count in the label to get the average
        % intensity of each cell of each channel
        averageTemplate = sumTemplate/pixels;
        averageTarget = sumTarget/pixels;
        averageThreshold = sumThreshold/pixels;
        
        % Sort cells into in and out of brightness threshold range
        if (averageTemplate > minThreshold && averageTemplate < maxThreshold && averageThreshold > minThreshold && averageThreshold < maxThreshold)
            cells = [cells; [averageTemplate averageTarget averageThreshold]];
            newregions(im_labels==i) = 1;
        else
            oldregions(im_labels==i) = 1;
        end
    end
    
    % Get new labels of cells based on the sorted cells
    [newlabels, newcount] = bwlabel(newregions);
    oldlabels = bwlabel(oldregions);
    
    % Create new boundaries from reduced label set
    newboundaries = bwboundaries(newlabels);
    oldboundaries = bwboundaries(oldlabels);
    
    % Show the channels as colors, and generate a heatmap of pixel-based
    % FRET signal
    if(showFigures && ~isempty(cells))
        
        % Creates color images
        zeroImage = zeros(512, 512);
        donorImage = cat(3, zeroImage, zeroImage, templateImage);
        acceptorImage = cat(3, zeroImage, thresholdImage, zeroImage);
        fretImage = cat(3, targetImage, zeroImage, zeroImage);
        heatmap = createHeatmap(newlabels, templateImage, targetImage, newcount);
        
        handle = figure('units','normalized','outerposition',[0 0 1 1]);
        
        plot1 = subplot_tight(1, 4, 1, [0.005]);
        imshow(donorImage);
        colorbar('Visible', 'off')
        title('Donor channel')

        plot2 = subplot_tight(1, 4, 2, [0.005]);
        imshow(acceptorImage);
        %showlabelinfo(newlabels, cells(:,2), newboundaries, oldboundaries, showText)
        colorbar('Visible', 'off')
        title('Acceptor channel')
        
        plot3 = subplot_tight(1, 4, 3, [0.005]);
        imshow(fretImage);
        %showlabelinfo(newlabels, cells(:,3), newboundaries, oldboundaries, showText)
        colorbar('Visible', 'off')
        title('FRET channel')
        
        plot4 = subplot_tight(1, 4, 4, [0.005]);
        imagesc(heatmap, [0,0.5]);
        axis image
        axis off
        colormap(hot)
        colorbar
        title('FRET heatmap')
        
        linkaxes([plot1, plot2, plot3, plot4],'xy') % Link all zoom and pan operations on the subplots
        
        if(showLabels)
            showlabelinfo(newlabels, cells(:,1), newboundaries, oldboundaries, plot1)
            showlabelinfo(newlabels, cells(:,1), newboundaries, oldboundaries, plot2)
            showlabelinfo(newlabels, cells(:,1), newboundaries, oldboundaries, plot3)
            showlabelinfo(newlabels, cells(:,1), newboundaries, oldboundaries, plot4)
        end
        
        addScalebars(handle, scale); % Add scalebars to all subplots
        
        % Add event listeners for zoom and pan operations on the figure, to
        % redraw the scalebars
        set(zoom(handle),'ActionPostCallback',@(obj, event) addScalebars(obj, scale));
        set(zoom(handle),'ActionPreCallback',@(obj, event) removeScalebars(obj));
        set(pan(handle),'ActionPostCallback',@(obj, event) addScalebars(obj, scale));
        set(pan(handle),'ActionPreCallback',@(obj, event) removeScalebars(obj));
        
        suptitle(filename)
        
        saveas(gcf,strcat('results/', filename, '_analysis.png'));
    end
    
end

function removeScalebars(handle)
    scaleTexts = findall(handle, 'tag', 'scaleText');
    scaleLines = findall(handle, 'tag', 'scaleLine');
    delete(scaleTexts)
    delete(scaleLines)
end

function addScalebars(handle, scale)
    % Loop through all axes on figure and add a scalebar
    axes = findall(handle,'type','axes');
    for i = 1:numel(axes)
        scalebar(axes(i), scale);
    end
end

% Creates a scalebar on axes based on the pixel/unit scale
function [line, scaleText] = scalebar(ax, scale)
    hold(ax, 'on')
    xlim = ax.XLim;
    ylim = ax.YLim;
    distance = xlim(2)-xlim(1);
    
    % Change scalebar size based on zoom level
    if(distance < 50)
        units = 50; % 50 um scalebar
    elseif(distance < 100)
        units = 100; % 100 um scalebar
    else
        units = 200; % 200 um scalebar
    end
    
    % Calculate the required scalebar size and position based on axes size and scale
    width = scale*units;
    x = xlim(1)+0.02*(xlim(2)-xlim(1));
    y = ylim(1)+0.02*(ylim(2)-ylim(1));
    offset = 0.02*(ylim(2)-ylim(1));
    
    % Show the scale line and text
    line = plot(ax, [x,x+width],[y,y],'Color','w','LineWidth',1, 'Tag', 'scaleLine');
    scaleText = text(ax, x+width/2,y+offset,strcat(num2str(units), '\mum'), 'Color', 'w', 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'Clipping', 'on', 'Tag', 'scaleText');
end

% Generates a pixel-based heatmap of the FRET signal
function [heatmap] = createHeatmap(labels, templateImage, targetImage, count)

    heatmap = zeros(512, 512);

    % Loops the cell labels
    for i = 1:count
        [r, c] = find(labels==i);
        rc = [r c];
        pixels = size(rc, 1);
        sumTemplate = 0;
        sumTarget = 0;
        fret = 0;

        % Loops the pixels of each label
        for j = 1:pixels
            valueTemplate = double(templateImage(rc(j,1), rc(j,2)));
            valueTarget = double(targetImage(rc(j,1), rc(j,2)));

            sumTemplate = sumTemplate + valueTemplate;
            sumTarget = sumTarget + valueTarget;

            % Get the FRET value and save it in the heatmap
            fret = valueTarget/(valueTemplate+valueTarget);
            heatmap(rc(j,1), rc(j,2)) = fret;
        end
    end
end

% Marks each label with a label number and average brightness
function showlabelinfo(newlabels, cells, newboundaries, oldboundaries, ax)
    global showText

    s = regionprops(newlabels, 'Centroid');
    
    hold(ax, 'on')
    visboundaries(ax, newboundaries, 'LineWidth', 0.1, 'EnhanceVisibility', false, 'Color', 'g');
    visboundaries(ax, oldboundaries, 'LineWidth', 0.1, 'EnhanceVisibility', false, 'Color', 'r');
    
    % Displays average cell intensity and label number
    if(showText)
        for k = 1:numel(s)
            c = s(k).Centroid;
            text(ax, c(1), c(2), sprintf('%d', int64(cells(k))), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', 'g', 'FontSize', 6, 'Clipping', 'on');
            text(ax, c(1), c(2)+5, sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', 'c', 'FontSize', 6, 'Clipping', 'on');
        end
    end
end