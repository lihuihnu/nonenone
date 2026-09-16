%% PLOT_GEOLOGY_MRST Visualize the proposed 20 x 20 x 5 geology with MRST.
%
% The script is intentionally independent from the reservoir executable.  It
% records a deliberately simple deterministic geology for the native
% PR/SW/CPA flow comparison and exports three separate, single-axes figures
% in English and Chinese.  Legends remain English in both variants.

clearvars;
close all;
clc;

mrstRoot = getenv('MRST_ROOT');
if isempty(mrstRoot)
    mrstRoot = 'D:\DOCUMENTS\研究生\MRST\mrst-2024b';
end
startupFile = fullfile(mrstRoot, 'startup.m');
assert(isfile(startupFile), 'MRST startup file was not found: %s', startupFile);
run(startupFile);

scriptDirectory = fileparts(mfilename('fullpath'));
caseDirectory = fileparts(scriptDirectory);
outputRoot = fullfile(caseDirectory, 'results', 'geology');
languages = {'en', 'zh'};
for languageIndex = 1:numel(languages)
    outputDirectory = fullfile(outputRoot, languages{languageIndex});
    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end
end

%% Grid and deterministic facies
cartDims = [20, 20, 5];
physicalDims = [1000, 600, 50]; % m
G = computeGeometry(cartGrid(cartDims, physicalDims));

[ii, jj, kk] = ind2sub(cartDims, (1:G.cells.num)');

% Facies codes:
% 1 background sand, 2 straight high-permeability channel,
% 3 shale baffle, 4 single baffle window.
facies = ones(G.cells.num, 1);

for cellIndex = 1:G.cells.num
    i = ii(cellIndex);
    j = jj(cellIndex);
    k = kk(cellIndex);

    if k == 3
        facies(cellIndex) = 3;
        inFlowWindow = i >= 9 && i <= 12 && j >= 9 && j <= 12;
        if inFlowWindow
            facies(cellIndex) = 4;
        end
        continue;
    end

    channelCenter = round(5 + 11 * (i - 2) / 17);
    if abs(j - channelCenter) <= 1
        facies(cellIndex) = 2;
    end
end

faciesCounts = accumarray(facies, 1, [4, 1]);
assert(G.cells.num == prod(cartDims), 'Unexpected MRST cell count.');
assert(all(faciesCounts > 0), 'At least one designed facies is absent.');
assert(faciesCounts(4) == 16, 'The baffle window must contain 16 cells.');

faciesNames = {
    'Background sand';
    'High-perm channel';
    'Shale baffle';
    'Flow window'};

% Color-blind-conscious qualitative colors.  Geometry is also distinguished
% by position, transparency, and well line style rather than color alone.
faciesColors = [
    0.82, 0.84, 0.86;  % background
    0.00, 0.45, 0.70;  % channel
    0.24, 0.25, 0.28;  % baffle
    0.00, 0.60, 0.44]; % window

permeabilityByFaciesMd = [
     80,  50,  5.00;
    250, 160, 18.00;
      2,   1,  0.05;
    150, 100, 10.00];
porosityByFacies = [0.19; 0.24; 0.10; 0.22];

kxMd = permeabilityByFaciesMd(facies, 1);
kyMd = permeabilityByFaciesMd(facies, 2);
kzMd = permeabilityByFaciesMd(facies, 3);
porosity = porosityByFacies(facies);

injectorIJ = [2, 5];
producerIJ = [19, 16];
injectorCells = sub2ind(cartDims, ...
    injectorIJ(1) * ones(2, 1), injectorIJ(2) * ones(2, 1), (1:2)');
producerCells = sub2ind(cartDims, ...
    producerIJ(1) * ones(2, 1), producerIJ(2) * ones(2, 1), (4:5)');

%% Accessible cell-wise source data
cellTable = table((1:G.cells.num)', ii, jj, kk, ...
    G.cells.centroids(:, 1), G.cells.centroids(:, 2), G.cells.centroids(:, 3), ...
    facies, string(faciesNames(facies)), kxMd, kyMd, kzMd, porosity, ...
    'VariableNames', {'cell', 'i', 'j', 'k', 'x_m', 'y_m', 'z_m', ...
    'facies_code', 'facies_name', 'kx_mD', 'ky_mD', 'kz_mD', 'porosity'});
writetable(cellTable, fullfile(outputRoot, 'geology_cells.csv'));

%% Export the three views in both language variants
for languageIndex = 1:numel(languages)
    language = languages{languageIndex};
    outputDirectory = fullfile(outputRoot, language);

    plotThreeDimensionalOverview(G, facies, faciesNames, faciesColors, ...
        injectorCells, producerCells, outputDirectory);
    plotReservoirChannelMap(G, cartDims, facies, faciesNames, faciesColors, ...
        injectorIJ, producerIJ, outputDirectory);
    plotBaffleLayerMap(G, cartDims, facies, faciesNames, faciesColors, ...
        injectorIJ, producerIJ, outputDirectory);
end

metadata = struct();
metadata.grid_cells = cartDims;
metadata.physical_dimensions_m = physicalDims;
metadata.cell_dimensions_m = physicalDims ./ cartDims;
metadata.facies_names = faciesNames;
metadata.permeability_mD = permeabilityByFaciesMd;
metadata.porosity = porosityByFacies;
metadata.injector_ij = injectorIJ;
metadata.injector_layers = [1, 2];
metadata.producer_ij = producerIJ;
metadata.producer_layers = [4, 5];
metadata.mrst_root = mrstRoot;
metadata.transformations = {
    'Deterministic facies assignment from logical cell indices';
    'No interpolation, smoothing, or stochastic realization';
    'One straight three-cell-wide high-permeability channel';
    'One 4-by-4 flow window in the middle shale layer';
    'Background cells rendered transparently only in the 3-D presentation';
    'The 3-D presentation uses five-fold vertical exaggeration';
    'Coordinate axes and internal cell-edge grids are hidden for presentation';
    'Only the outer geological-grid boundary frame is retained'};
metadataFile = fopen(fullfile(outputRoot, 'geology_manifest.json'), 'w');
assert(metadataFile >= 0, 'Could not create geology manifest.');
cleanup = onCleanup(@() fclose(metadataFile));
fprintf(metadataFile, '%s\n', jsonencode(metadata, PrettyPrint=true));

fprintf('MRST geology figures written to: %s\n', outputRoot);

%% Local plotting functions
function plotThreeDimensionalOverview(G, facies, faciesNames, colors, ...
    injectorCells, producerCells, outputDirectory)
    fig = newFigure([17.0, 12.0]);
    ax = axes(fig);
    hold(ax, 'on');

    plotGrid(G, facies == 1, 'FaceColor', colors(1, :), ...
        'FaceAlpha', 0.035, 'EdgeColor', 'none');
    plotGrid(G, facies == 2, 'FaceColor', colors(2, :), ...
        'FaceAlpha', 0.88, 'EdgeColor', 'none');
    plotGrid(G, facies == 3, 'FaceColor', colors(3, :), ...
        'FaceAlpha', 0.34, 'EdgeColor', 'none');
    plotGrid(G, facies == 4, 'FaceColor', colors(4, :), ...
        'FaceAlpha', 0.96, 'EdgeColor', [0.00, 0.30, 0.22], ...
        'LineWidth', 0.45);
    drawThreeDimensionalBoundary(ax, G);

    injector = G.cells.centroids(injectorCells, :);
    producer = G.cells.centroids(producerCells, :);
    injectorStem = [injector(1, 1:2), -12; injector(1, 1:2), 20];
    producerStem = [producer(1, 1:2), -12; producer(1, 1:2), 50];
    plot3(ax, injectorStem(:, 1), injectorStem(:, 2), injectorStem(:, 3), '-', ...
        'Color', [0.50, 0.12, 0.72], 'LineWidth', 5.0);
    plot3(ax, injector(:, 1), injector(:, 2), injector(:, 3), 'o', ...
        'MarkerFaceColor', [0.50, 0.12, 0.72], 'MarkerEdgeColor', 'white', ...
        'MarkerSize', 7.5, 'LineWidth', 1.1);
    plot3(ax, producerStem(:, 1), producerStem(:, 2), producerStem(:, 3), '--', ...
        'Color', [0.05, 0.05, 0.05], 'LineWidth', 4.5);
    plot3(ax, producer(:, 1), producer(:, 2), producer(:, 3), 's', ...
        'MarkerFaceColor', 'white', 'MarkerEdgeColor', [0.05, 0.05, 0.05], ...
        'MarkerSize', 7.2, 'LineWidth', 1.2);

    axis(ax, 'equal');
    xlim(ax, [0, 1000]);
    ylim(ax, [0, 600]);
    zlim(ax, [-14, 50]);
    view(ax, 42, 27);
    set(ax, 'ZDir', 'reverse');
    daspect(ax, [1, 1, 0.2]);
    applyImageOnlyStyle(ax);
    addFaciesLegend(ax, faciesNames, colors, true);
    exportFigurePair(fig, outputDirectory, '01_geology_3d');
end

function plotReservoirChannelMap(G, cartDims, facies, faciesNames, colors, ...
    injectorIJ, producerIJ, outputDirectory)
    fig = newFigure([15.0, 10.8]);
    ax = axes(fig);
    hold(ax, 'on');

    layerCells = find(repelem((1:cartDims(3))', prod(cartDims(1:2))) == 1);
    for code = 1:2
        cells = layerCells(facies(layerCells) == code);
        plotGrid(G, cells, 'FaceColor', colors(code, :), ...
            'EdgeColor', 'none');
    end
    axis(ax, 'equal');
    axis(ax, 'tight');
    view(ax, 0, 90);
    set(ax, 'SortMethod', 'childorder');
    drawPlanBoundary(ax, G);
    markProjectedWells(ax, G, cartDims, injectorIJ, producerIJ);
    applyImageOnlyStyle(ax);
    addFaciesLegend(ax, faciesNames(1:2), colors(1:2, :), true);
    exportFigurePair(fig, outputDirectory, '02_channel_map');
end

function plotBaffleLayerMap(G, cartDims, facies, faciesNames, colors, ...
    injectorIJ, producerIJ, outputDirectory)
    fig = newFigure([15.0, 10.8]);
    ax = axes(fig);
    hold(ax, 'on');

    layerStart = prod(cartDims(1:2)) * 2 + 1;
    layerCells = (layerStart:(layerStart + prod(cartDims(1:2)) - 1))';
    for code = 3:4
        cells = layerCells(facies(layerCells) == code);
        plotGrid(G, cells, 'FaceColor', colors(code, :), ...
            'EdgeColor', 'none');
    end
    axis(ax, 'equal');
    axis(ax, 'tight');
    view(ax, 0, 90);
    set(ax, 'SortMethod', 'childorder');
    drawPlanBoundary(ax, G);
    markProjectedWells(ax, G, cartDims, injectorIJ, producerIJ);
    applyImageOnlyStyle(ax);
    addFaciesLegend(ax, faciesNames(3:4), colors(3:4, :), true);
    exportFigurePair(fig, outputDirectory, '03_baffle_layer');
end

function markProjectedWells(ax, G, cartDims, injectorIJ, producerIJ)
    injectorCell = sub2ind(cartDims, injectorIJ(1), injectorIJ(2), 1);
    producerCell = sub2ind(cartDims, producerIJ(1), producerIJ(2), 1);
    injectorXY = G.cells.centroids(injectorCell, 1:2);
    producerXY = G.cells.centroids(producerCell, 1:2);
    scatter(ax, injectorXY(1), injectorXY(2), 92, ...
        'MarkerFaceColor', [0.50, 0.12, 0.72], 'MarkerEdgeColor', 'white', ...
        'LineWidth', 1.2, 'Clipping', 'off');
    scatter(ax, producerXY(1), producerXY(2), 82, 's', ...
        'MarkerFaceColor', 'white', 'MarkerEdgeColor', [0.05, 0.05, 0.05], ...
        'LineWidth', 1.4, 'Clipping', 'off');
end

function drawThreeDimensionalBoundary(ax, G)
    lower = min(G.nodes.coords, [], 1);
    upper = max(G.nodes.coords, [], 1);
    corners = [
        lower(1), lower(2), lower(3);
        upper(1), lower(2), lower(3);
        upper(1), upper(2), lower(3);
        lower(1), upper(2), lower(3);
        lower(1), lower(2), upper(3);
        upper(1), lower(2), upper(3);
        upper(1), upper(2), upper(3);
        lower(1), upper(2), upper(3)];
    edges = [
        1, 2; 2, 3; 3, 4; 4, 1;
        5, 6; 6, 7; 7, 8; 8, 5;
        1, 5; 2, 6; 3, 7; 4, 8];
    for edgeIndex = 1:size(edges, 1)
        nodes = corners(edges(edgeIndex, :), :);
        plot3(ax, nodes(:, 1), nodes(:, 2), nodes(:, 3), '-', ...
            'Color', [0.26, 0.27, 0.29], 'LineWidth', 0.78);
    end
end

function drawPlanBoundary(ax, G)
    lower = min(G.nodes.coords(:, 1:2), [], 1);
    upper = max(G.nodes.coords(:, 1:2), [], 1);
    x = [lower(1), upper(1), upper(1), lower(1), lower(1)];
    y = [lower(2), lower(2), upper(2), upper(2), lower(2)];
    z = max(zlim(ax)) * ones(size(x));
    plot3(ax, x, y, z, '-', 'Color', [0.26, 0.27, 0.29], ...
        'LineWidth', 0.92, 'Clipping', 'off');
end

function addFaciesLegend(ax, names, colors, includeWells)
    hold(ax, 'on');
    handles = gobjects(numel(names) + 2 * includeWells, 1);
    labels = cell(numel(handles), 1);
    for index = 1:numel(names)
        handles(index) = plot3(ax, nan, nan, nan, 's', ...
            'MarkerFaceColor', colors(index, :), ...
            'MarkerEdgeColor', [0.18, 0.18, 0.18], ...
            'MarkerSize', 7.2, 'LineStyle', 'none');
        labels{index} = names{index};
    end
    if includeWells
        offset = numel(names);
        handles(offset + 1) = plot3(ax, nan, nan, nan, 'o', ...
            'MarkerFaceColor', [0.50, 0.12, 0.72], ...
            'MarkerEdgeColor', 'white', 'MarkerSize', 7.2, ...
            'LineStyle', 'none');
        handles(offset + 2) = plot3(ax, nan, nan, nan, 's', ...
            'MarkerFaceColor', 'white', ...
            'MarkerEdgeColor', [0.05, 0.05, 0.05], 'MarkerSize', 7.2, ...
            'LineStyle', 'none');
        labels{offset + 1} = 'CO_2 injector';
        labels{offset + 2} = 'Producer';
    end
    lgd = legend(ax, handles, labels, 'Location', 'northwest', ...
        'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 7.5, 'Box', 'on');
    lgd.ItemTokenSize = [14, 10];
end

function fig = newFigure(sizeCm)
    fig = figure('Color', 'white', 'Units', 'centimeters', ...
        'Position', [2, 2, sizeCm(1), sizeCm(2)], ...
        'Visible', 'off');
end

function applyImageOnlyStyle(ax)
    grid(ax, 'off');
    set(ax, 'Visible', 'off', 'Box', 'off', ...
        'Position', [0.015, 0.015, 0.97, 0.97]);
end

function exportFigurePair(fig, outputDirectory, stem)
    pngPath = fullfile(outputDirectory, [stem, '.png']);
    pdfPath = fullfile(outputDirectory, [stem, '.pdf']);
    exportgraphics(fig, pngPath, 'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    close(fig);
end
