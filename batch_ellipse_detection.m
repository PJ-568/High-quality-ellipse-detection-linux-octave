% Batch Ellipse Detection Script/Function
% This script processes all images in an input folder, detects ellipses,
% and saves results to an output folder.
%
% Usage as script:
% 1. Modify the parameters below (input_folder, output_folder, Tac, Tr, specified_polarity).
% 2. Run this script in Octave: octave batch_ellipse_detection.m
%
% Usage as function from Octave command line:
%   batch_ellipse_detection(input_folder, output_folder, Tac, Tr, specified_polarity)
%
% Usage from bash (via octave --eval):
%   octave --eval "batch_ellipse_detection('pics', 'results', 165, 0.6, 0)"
%
% Results:
% - For each image, a result image with ellipses drawn is saved as [image_name]_result.jpg
% - A text file with ellipse parameters is saved as [image_name]_ellipses.txt
% - A summary CSV file is generated with statistics for all processed images.

function [] = batch_ellipse_detection(input_folder, output_folder, Tac, Tr, specified_polarity)

% If called with no arguments, use default values
if nargin < 1
    input_folder = 'pics';
end
if nargin < 2
    output_folder = 'results';
end
if nargin < 3
    Tac = 165;          % elliptic angular coverage (completeness degree)
end
if nargin < 4
    Tr = 0.6;           % ratio of support inliers
end
if nargin < 5
    specified_polarity = 0;  % 0: all ellipses, 1: positive polarity, -1: negative polarity
end

% Image file extensions to process
image_extensions = {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff'};

% Create output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Get list of image files
image_files = {};
for ext = image_extensions
    files = dir(fullfile(input_folder, ['*' ext{1}]));
    for i = 1:length(files)
        image_files{end+1} = fullfile(input_folder, files(i).name);
    end
end

% Sort alphabetically
image_files = sort(image_files);

if isempty(image_files)
    error('No image files found in %s with extensions: %s', input_folder, strjoin(image_extensions, ', '));
end

fprintf('Found %d image(s) to process.\n', length(image_files));

%% Initialize summary data
summary_data = cell(length(image_files)+1, 6);
summary_data{1,1} = 'Image';
summary_data{1,2} = 'Width';
summary_data{1,3} = 'Height';
summary_data{1,4} = 'Num_Ellipses';
summary_data{1,5} = 'Processing_Time(s)';
summary_data{1,6} = 'Result_File';

%% Process each image
for img_idx = 1:length(image_files)
    image_path = image_files{img_idx};
    [~, image_name, image_ext] = fileparts(image_path);
    
    fprintf('\nProcessing image %d/%d: %s\n', img_idx, length(image_files), image_name);
    
    % Read image
    try
        I = imread(image_path);
    catch
        fprintf('Warning: Could not read %s, skipping.\n', image_path);
        continue;
    end
    
    % Get image dimensions
    [height, width, ~] = size(I);
    
    % Start timer
    t_start = tic;
    
    % Detect ellipses
    try
        [ellipses, ~, ~] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);
        num_ellipses = size(ellipses, 1);
    catch err
        fprintf('Error during ellipse detection: %s\n', err.message);
        num_ellipses = 0;
        ellipses = [];
    end
    
    % Stop timer
    processing_time = toc(t_start);
    
    % Convert angle from radians to degrees for display
    if ~isempty(ellipses)
        ellipses(:,5) = ellipses(:,5) ./ pi * 180;
    end
    
    % Save ellipse parameters to text file
    ellipses_filename = fullfile(output_folder, [image_name '_ellipses.txt']);
    fid = fopen(ellipses_filename, 'w');
    if fid ~= -1
        fprintf(fid, 'Ellipse parameters for %s\n', image_name);
        fprintf(fid, 'Format: center_x, center_y, a (semi-major), b (semi-minor), phi (degrees)\n');
        fprintf(fid, 'Number of ellipses: %d\n\n', num_ellipses);
        for i = 1:size(ellipses, 1)
            fprintf(fid, '%.2f, %.2f, %.2f, %.2f, %.2f\n', ...
                ellipses(i,1), ellipses(i,2), ellipses(i,3), ellipses(i,4), ellipses(i,5));
        end
        fclose(fid);
    else
        fprintf('Warning: Could not write ellipse data to %s\n', ellipses_filename);
    end
    
    % Draw ellipses on image and save result
    if ~isempty(ellipses)
        % Create figure without displaying it
        fig = figure('Visible', 'off');
        imshow(I);
        hold on;
        
        % Draw each ellipse
        th = 0:pi/180:2*pi;
        for i = 1:size(ellipses, 1)
            Semi_major = ellipses(i,3);
            Semi_minor = ellipses(i,4);
            x0 = ellipses(i,1);
            y0 = ellipses(i,2);
            Phi = ellipses(i,5) * pi / 180; % Convert back to radians
            
            x = x0 + Semi_major*cos(Phi)*cos(th) - Semi_minor*sin(Phi)*sin(th);
            y = y0 + Semi_minor*cos(Phi)*sin(th) + Semi_major*sin(Phi)*cos(th);
            
            plot(x, y, 'r', 'LineWidth', 2);
        end
        
        % Configure axes
        axis on;
        set(gca, 'XTick', [], 'YTick', []);
        axis ij;
        axis equal;
        axis([0 width 0 height]);
        
        % Save figure
        result_filename = fullfile(output_folder, [image_name '_result.jpg']);
        print(fig, result_filename, '-djpeg', '-r150');
        close(fig);
    else
        % If no ellipses detected, just save the original image as result
        result_filename = fullfile(output_folder, [image_name '_result.jpg']);
        imwrite(I, result_filename);
    end
    
    % Update summary data
    summary_data{img_idx+1,1} = image_name;
    summary_data{img_idx+1,2} = width;
    summary_data{img_idx+1,3} = height;
    summary_data{img_idx+1,4} = num_ellipses;
    summary_data{img_idx+1,5} = processing_time;
    summary_data{img_idx+1,6} = [image_name '_result.jpg'];
    
    fprintf('  Detected %d ellipse(s) in %.3f seconds\n', num_ellipses, processing_time);
end

%% Save summary CSV file
summary_filename = fullfile(output_folder, 'processing_summary.csv');
fid = fopen(summary_filename, 'w');
if fid ~= -1
    % Write header
    fprintf(fid, '%s,%s,%s,%s,%s,%s\n', ...
        'Image', 'Width', 'Height', 'Num_Ellipses', 'Processing_Time(s)', 'Result_File');
    
    % Write data rows
    for i = 2:size(summary_data, 1)
        if ~isempty(summary_data{i,1})
            fprintf(fid, '%s,%d,%d,%d,%.3f,%s\n', ...
                summary_data{i,1}, summary_data{i,2}, summary_data{i,3}, ...
                summary_data{i,4}, summary_data{i,5}, summary_data{i,6});
        end
    end
    fclose(fid);
    fprintf('\nSummary saved to: %s\n', summary_filename);
else
    fprintf('Warning: Could not write summary CSV file\n');
end

fprintf('\nBatch processing completed!\n');
fprintf('Results saved in: %s\n', output_folder);

end