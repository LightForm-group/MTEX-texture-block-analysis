%% Specify Crystal and Specimen Symmetries

% crystal symmetry with alpha first
CS = {... 
   'notIndexed',...
   crystalSymmetry('6/mmm', [2.954 2.954 4.729], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', 'light blue'),...
   crystalSymmetry('m-3m', [3.192 3.192 3.192], 'mineral', 'Titanium cubic', 'color', 'light green')};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names

% path to files
pname = '/Users/mbcx9cd4/Documents/MATLAB/ebsd/MTEX-texture-block-analysis/';

% which files to be imported
sample_name = 'Ti64_915C_87pct_Sample1_TD-ND_Plane_Beta'
data_path = strcat('Data/',sample_name,'/',sample_name,'_Data.crc')

analysis_path = strcat('Analysis/',sample_name,'/') % path for saving the data

% which files to be imported
fname = [pname data_path]

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','crc',...
'convertEuler2SpatialReferenceFrame')

%% Reduce the Data

ebsd = reduce(ebsd)   % take every second pixel horiz. and vert.

%% Rotate the data

rot = rotation('Euler', 90*degree, 90*degree, 0*degree);

ebsd = rotate(ebsd,rot,'keepXY'); % rotate the orientation data
ebsd = rotate(ebsd,88.5*degree,'keepEuler') % rotate the spatial data

fprintf('Note, the (x,y) origin on the map will have changed and x or y could be negative!')

%% Set whether figures are visible

visible = 'off'

%% Plot the IPF colour map

phase = 'beta'
outputFileName = strcat(analysis_path,sample_name,'_IPF_map_entire_region')
IPF_map_plot(phase, ebsd, outputFileName, visible)

%% Plot the pole figures for the whole sample

phase = 'beta'
ori = ebsd('Titanium cubic').orientations
contour_step = 0.25
pf_max = 3.0
outputFileName = strcat(analysis_path,sample_name,'_PF_entire_region')
pole_figure_plot(phase, ori, CS, contour_step, pf_max, outputFileName, visible);

%% Change the specimen symmetry for calculating the ODF

ori.SS = specimenSymmetry('orthorhombic')

%% Calculate an automatic halfwidth for the ODF (not needed in this case)

% if the orientations are spatially independant...
% psi=calcKernel(ebsd('Ti-Hex').orientations);

% if the EBSD measurements are not rough texture measurements and are spatially dependant (more than one point per grain), then perform grain reconstruction and estimate the halfwidth from the grains...
% grain reconstruction (default is 10 degrees so could just use calcGrains(ebsd))...
% grains = calcGrains(ebsd,'angle',10*degree);
% correct for small grains...
% grains = grains(grains.grainSize>5);
% compute optimal halfwidth from the meanorientations of grains...
% psi = calcKernel(grains('Zirc-alloy4').meanOrientation)

% HALF_WIDTH = psi

%% Calculate the ODF using the optimised halfwidth

%odf = calcDensity(ori,'kernel',psi);
odf = calcDensity(ori,'halfwidth',5*degree);

%% Calculate the Texture Index for the whole sample

TEXTURE_INDEX = textureindex(odf)
MAX = max(odf)

%% Plot the ODF slices without contouring for the whole sample

odf_max = 5.0
outputFileName = strcat(analysis_path,sample_name,'_ODF_entire_region')
specSym = 'orthorhombic'
ODF_plot(phase, odf, odf_max, outputFileName, specSym, visible)

%% Change the specimen symmetry back to triclinic

ori.SS = specimenSymmetry('triclinic')

%% Crop the map into a rectangle (avoiding the unindexed edges of the sample)
% note, in this case x is vertical starting at 0 and runs downwards as
% negative values, y is horizontal starting at 0 and runs left-to-right as
% positive values. To show the map comment out close(IPF_map) in 
% IPF_map_plot.m then choose the (x,y) points by clicking on the figure.

x_top = -130
y_left = 220
x_bottom = -4640
y_right = 20920

x_width = x_bottom-x_top
y_width = y_right-y_left

region = [x_top, y_left, x_width, y_width]; % note, region is defined as x,y origin and an x,y width which is added onto the origin
condition = inpolygon(ebsd,region); % points located within region
ebsd_cropped = ebsd(condition);
ori_cropped = ebsd_cropped('Titanium cubic').orientations

% change the specimen symmetry to triclinic for plotting the PFs
ori_cropped.SS = specimenSymmetry('triclinic')

% plot the IPF map for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_IPF_map_cropped')
IPF_map_plot(phase, ebsd_cropped, outputFileName, visible)

% plot the pole figures for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_PF_cropped')
pole_figure_plot(phase, ori_cropped, CS, contour_step, pf_max, outputFileName, visible);

% change the specimen symmetry to orthorhombic for plotting the ODF
ori_cropped.SS = specimenSymmetry('orthorhombic')

% calculate an automatic halfwidth for the ODF for the cropped compression sample
% psi=calcKernel(ori_cropped);
% HALF_WIDTH = psi

% calculate the ODF using the optimised halfwidth for the cropped compression sample
% odf_cropped = calcDensity(ori_cropped,'kernel',psi);
odf_cropped = calcDensity(ori_cropped,'halfwidth',5*degree);

% calculate the Texture Index for the cropped compression sample
TEXTURE_INDEX = textureindex(odf_cropped)
MAX = max(odf_cropped)

% plot the ODF slices without contouring for the cropped compression sample
outputFileName = strcat(analysis_path,sample_name,'_ODF_cropped')
specSym = 'orthorhombic'
ODF_plot(phase, odf_cropped, odf_max, outputFileName, specSym, visible)

%% Change the specimen symmetry back to triclinic

ori_cropped.SS = specimenSymmetry('triclinic')

%% Choose whether to slice the full map or the cropped map here

% use for cropped map
ebsd = ebsd_cropped
x_origin = x_top
y_origin = y_left

% use for entire map
% ebsd = ebsd
% x_origin = 0
% y_origin = 0

%% Slice the entire map or the cropped map

num_squares_x = 9; % number of squares to cut the map into in x (vertical)
num_squares_y = 43; % number of squares to cut the map into in y (horizontal)

% define the size of the EBSD map
ebsd_grid = ebsd.gridify;
ebsd_shape = size(ebsd_grid.id);
original_y = ebsd_shape(1);
original_x = ebsd_shape(2);
stepSize = ebsd_grid.dx;

x_min = (sqrt(x_origin * x_origin)/stepSize);
x_max = original_x + (sqrt(x_origin * x_origin)/stepSize);
x_length = x_max - x_min;

y_min = (sqrt(y_origin * y_origin)/stepSize);
y_max = original_y + (sqrt(y_origin * y_origin)/stepSize);
y_length = y_max - y_min;

% for splitting into squares along x
x_width = floor(x_length / num_squares_x); % round to nearest integer
x_axis = (1:num_squares_x);

% for splitting into squares along y
y_width = floor(y_length / num_squares_y); % round to nearest integer
y_axis = (1:num_squares_y);

cutmap = containers.Map('KeyType', 'int32', 'ValueType', 'any'); % creates an empty Map object

square_number = 1

for strip_index_x = 0:num_squares_x-1
    
    for strip_index_y = 0:num_squares_y-1
        % separate the map section

        % set out the coordinates for the edge of the region
        % note, region is defined as x,y origin and an x,y width which is added onto the origin

        % if splitting into squares and x is negative
        x_min_square = strip_index_x * x_width + x_min;
        y_min_square = strip_index_y * y_width + y_min;
        region = [-x_min_square*stepSize, y_min_square*stepSize, -x_width*stepSize, y_width*stepSize]
        
        % if splitting into squares and x is positive
        %x_min_square = strip_index_x * x_width + x_min;
        %y_min_square = strip_index_y * y_width + y_min;
        %region = [x_min_square*stepSize, y_min_square*stepSize, x_width*stepSize, y_length*stepSize]
        
        % Cut the EBSD map
        condition = inpolygon(ebsd,region); % points located within region
        ebsd_square = ebsd(condition); % create ebsd map for region
        cutmap(square_number) = ebsd_square; % store square in Map object with index
        ebsd_cutmap = cutmap(square_number); % read out ebsd_cutmap from the Map object
        
        % plot the IPF map to check the slices
        outputFileName = strcat(analysis_path,sample_name,'_IPF_map_square_',num2str(square_number))
        IPF_map_plot(phase, ebsd_cutmap, outputFileName, visible)  
      
        square_number = square_number + 1
    end    
end 

%% Analyse and plot the sliced data to see how the texture components change along the length

% define the crystal system for the texture components
cs = ebsd('Titanium cubic').CS;

% define the maximum possible misorientation
misorientation = 5

% Define a texture component for the hexagonal phase
% basal_ND = symmetrise(orientation.byMiller([0 0 0 1],[1 0 -1 0],cs),'unique'); % define component with directions
% basal_ND = symmetrise(orientation.byEuler(0*degree,0*degree,0*degree,cs),'unique') % define component with Euler angles
% tilted_ND = symmetrise(orientation.byEuler(0*degree,30*degree,0*degree,cs),'unique') % define component with Euler angles
basal_TD = symmetrise(orientation.byEuler(0*degree,90*degree,0*degree,cs),'unique') % define component with Euler angles
basal_RD = symmetrise(orientation.byEuler(90*degree,90*degree,0*degree,cs),'unique') % define component with Euler angles

% Define a texture component for the cubic phase
rotated_cube = orientation.byMiller([0 0 1],[0 1 1],cs); % define component with directions
% Define a texture fibre for the cubic phase
gamma_fibre = fibre(Miller(1,1,1,cs),xvector);
alpha_fibre = fibre(Miller(1,1,0,cs),yvector);

num_squares = num_squares_x * num_squares_y

for square_index = 1:num_squares
    
    ebsd_cutmap = cutmap(square_index); % read out ebsd_cutmap from the Map object
   
    % plot the IPF map, pole figures and odf slices
    outputFileName = strcat(analysis_path,sample_name,'_IPF_map_square_',num2str(square_index))
    IPF_map_plot(phase, ebsd_cutmap, outputFileName, visible)
    
    ori_square = ebsd_cutmap('Titanium cubic').orientations
    % change the specimen symmetry to triclinic for plotting the PFs
    ori_square.SS = specimenSymmetry('triclinic')
    outputFileName = strcat(analysis_path,sample_name,'_PF_square_',num2str(square_index))
    [maxval] = pole_figure_plot(phase, ori_square, CS, contour_step, pf_max, outputFileName, visible);
    PF_001_max(square_index) = maxval(1);
    PF_110_max(square_index) = maxval(2);
    PF_111_max(square_index) = maxval(3);
    
    outputFileName = strcat(analysis_path,sample_name,'_ODF_square_',num2str(square_index))
    % change the specimen symmetry to orthorhombic for plotting the ODF
    ori_square.SS = specimenSymmetry('orthorhombic')
    %psi=calcKernel(ori_square);
    %HALF_WIDTH = psi
    %odf_square = calcDensity(ori_square,'kernel',psi);
    odf_square = calcDensity(ori_square,'halfwidth',5*degree);
    ODF_plot(phase, odf_square, odf_max, outputFileName, specSym, visible)
     
    % calculate texture index
    TEXTURE_INDEX_square(square_index) = textureindex(odf_square)
    
    % change the specimen symmetry to triclinic for plotting the PFs
    ori_square.SS = specimenSymmetry('triclinic')
    
    % calculate strength of ODF maxima
    [odf_square_max(square_index),ori_square_max(square_index)]= max(odf_square)
    
    % calculate misorientation of ODF maxima wrt 0002 in CD
    %mori = (orientation.byEuler(0*degree,0*degree,0*degree,cs))*ori_square_max
    %misorientation_ODF_max = angle(mori)/degree

    % calculate PHI angle of ODF maxima
    phi1(square_index) = ori_square_max(square_index).phi1
    PHI(square_index) = ori_square_max(square_index).Phi
    phi2(square_index) = ori_square_max(square_index).phi2
    
    % seperate the texture components and calculate the volume fractions
    %total_volume = length(ebsd_cutmap) % calculate the total volume as the number of points in the map
    
    % seperate a texture component and calculate the volume fraction
    rotated_cube_volume_fraction(square_index) = volume(odf_square, rotated_cube, misorientation*degree);
    alpha_fibre_volume_fraction(square_index) = volume(odf_square, alpha_fibre, misorientation*degree);
    gamma_fibre_volume_fraction(square_index) = volume(odf_square, gamma_fibre, misorientation*degree);
end

%% Open and write to file to save the different texture strength values

fileTS = fopen(fullfile(analysis_path, strcat(sample_name,'_texture_strength_',num2str(num_squares),'.txt')),'w');
fprintf(fileTS, 'Square Index \t Texture Index \t ODF Max \t phi1 Angle of ODF Max \t PHI Angle of ODF Max \t phi2 Angle of ODF Max \t {001} PF Max \t {110} PF Max \t {111} PF Max \t Rotated Cube Volume Fraction \t Alpha Fibre Volume Fraction \t Gamma Fibre Volume Fraction \n');

for square_index = 1:num_squares    
    % write the texture strength values to file
    fprintf(fileTS, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', square_index, TEXTURE_INDEX_square(square_index), odf_square_max(square_index), rad2deg(phi1(square_index)), rad2deg(PHI(square_index)), rad2deg(phi2(square_index)), PF_001_max(square_index), PF_110_max(square_index), PF_111_max(square_index), rotated_cube_volume_fraction(square_index), alpha_fibre_volume_fraction(square_index), gamma_fibre_volume_fraction(square_index))
end

% close any open files
fclose(fileTS);

%% Plot the texture variation

square_index = 1:num_squares

TEXTURE_INDEX_figure = figure();
hold on
plot(square_index,TEXTURE_INDEX_square,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Texture Index')
hold off
legend('Texture Index')
saveas (TEXTURE_INDEX_figure, strcat(analysis_path,sample_name,'_TI_line_plot_',num2str(num_squares), '.png'));

odf_max_figure = figure();
hold on
plot(square_index,odf_square_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('ODF Maximum')
hold off
legend('ODF Maximum')
saveas (odf_max_figure, strcat(analysis_path,sample_name,'_odf_max_line_plot_',num2str(num_squares), '.png'));

phi1_figure = figure();
hold on
plot(square_index,rad2deg(phi1),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Phi1 Angle')
hold off
legend('Phi1 Angle')
saveas (phi1_figure, strcat(analysis_path,sample_name,'_phi1_line_plot_',num2str(num_squares), '.png'));

PHI_figure = figure();
hold on
plot(square_index,rad2deg(PHI),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('PHI Angle')
hold off
legend('PHI Angle')
saveas (PHI_figure, strcat(analysis_path,sample_name,'_PHI_line_plot_',num2str(num_squares), '.png'));

phi2_figure = figure();
hold on
plot(square_index,rad2deg(phi2),'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Phi2 Angle')
hold off
legend('Phi2 Angle')
saveas (phi2_figure, strcat(analysis_path,sample_name,'_phi2_line_plot_',num2str(num_squares), '.png'));

PF_001_figure = figure();
hold on
plot(square_index,PF_001_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{001} Pole Figure Max.')
hold off
legend('{001} Pole Figure Max.')
saveas (PF_001_figure, strcat(analysis_path,sample_name,'_PF_001_line_plot_',num2str(num_squares), '.png'));

PF_110_figure = figure();
hold on
plot(square_index,PF_110_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{110} Pole Figure Max.')
hold off
legend('{110} Pole Figure Max.')
saveas (PF_110_figure, strcat(analysis_path,sample_name,'_PF_110_line_plot_',num2str(num_squares), '.png'));

PF_111_figure = figure();
hold on
plot(square_index,PF_111_max,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('{111} Pole Figure Max.')
hold off
legend('{111} Pole Figure Max.')
saveas (PF_111_figure, strcat(analysis_path,sample_name,'_PF_111_line_plot_',num2str(num_squares), '.png'));

rotated_cube_volume_fraction_figure = figure();
hold on
plot(square_index,rotated_cube_volume_fraction*100,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('R.C. Volume Fraction (%)')
hold off
legend('R.C. Volume Fraction')
saveas (rotated_cube_volume_fraction_figure, strcat(analysis_path,sample_name,'_rotated_cube_volume_fraction_line_plot_',num2str(num_squares), '.png'));

alpha_fibre_volume_fraction_figure = figure();
hold on
plot(square_index,alpha_fibre_volume_fraction*100,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Alpha Fibre Volume Fraction (%)')
hold off
legend('Alpha Fibre Volume Fraction')
saveas (alpha_fibre_volume_fraction_figure, strcat(analysis_path,sample_name,'_alpha_fibre_volume_fraction_line_plot_',num2str(num_squares), '.png'));

gamma_fibre_volume_fraction_figure = figure();
hold on
plot(square_index,gamma_fibre_volume_fraction*100,'Color',[1,0,0],'lineWidth',2) % red;
xlabel('Slice Number')
ylabel('Gamma Fibre Volume Fraction (%)')
hold off
legend('Gamma Fibre Volume Fraction')
saveas (gamma_fibre_volume_fraction_figure, strcat(analysis_path,sample_name,'_gamma_fibre_volume_fraction_line_plot_',num2str(num_squares), '.png'));
