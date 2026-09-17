function [fig, outputFiles] = plot_geological_model(outputDirectory)
%PLOT_GEOLOGICAL_MODEL Plot the structured geometry of scw_kerogen_lmh_2d.
%   PLOT_GEOLOGICAL_MODEL() reads the authoritative grid dimensions from
%   ../scw_kerogen_common/benchmark_2d_common.hpp, draws the 60-by-20-by-1
%   homogeneous reservoir, marks the SCW injector and producer, and exports
%   publication-ready PNG, PDF, SVG, and MATLAB FIG files.
%
%   [FIG, FILES] = PLOT_GEOLOGICAL_MODEL(OUTPUTDIRECTORY) selects a custom
%   export directory and returns the figure handle and exported file paths.

    arguments
        outputDirectory (1, 1) string = ""
    end

    scriptDirectory = fileparts(mfilename('fullpath'));
    caseDirectory = fileparts(scriptDirectory);
    commonConfig = fullfile(caseDirectory, '..', ...
        'scw_kerogen_common', 'benchmark_2d_common.hpp');

    if ~isfile(commonConfig)
        error('Grid configuration not found: %s', commonConfig);
    end
    if strlength(outputDirectory) == 0
        outputDirectory = fullfile(caseDirectory, 'figures', ...
            'geological_model');
    end
    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    sourceText = fileread(commonConfig);
    gridDimensions = [ ...
        readNumericConstant(sourceText, 'nx'), ...
        readNumericConstant(sourceText, 'ny'), ...
        readNumericConstant(sourceText, 'nz')];
    physicalDimensions = [ ...
        readNumericConstant(sourceText, 'lx'), ...
        readNumericConstant(sourceText, 'ly'), ...
        readNumericConstant(sourceText, 'lz')];

    nx = gridDimensions(1);
    ny = gridDimensions(2);
    nz = gridDimensions(3);
    lx = physicalDimensions(1);
    ly = physicalDimensions(2);
    lz = physicalDimensions(3);

    if any(gridDimensions < 1) || any(mod(gridDimensions, 1) ~= 0)
        error('Grid dimensions must be positive integers.');
    end
    if any(physicalDimensions <= 0)
        error('Physical dimensions must be positive.');
    end

    % The production case defines primaryCentreJ = ny/2 - 1 (zero based).
    primaryCentreJ = ny / 2 - 1;
    injectorI = 0;
    producerI = nx - 1;
    cellSize = physicalDimensions ./ gridDimensions;
    injector = [(injectorI + 0.5) * cellSize(1), ...
        (primaryCentreJ + 0.5) * cellSize(2), 0.5 * cellSize(3)];
    producer = [(producerI + 0.5) * cellSize(1), ...
        (primaryCentreJ + 0.5) * cellSize(2), 0.5 * cellSize(3)];

    paleBlue = [0.78, 0.87, 0.94];
    gridBlue = [0.43, 0.59, 0.72];
    darkBlue = [0.12, 0.31, 0.48];
    injectorBlue = [0.16, 0.55, 0.72];
    producerBlue = [0.08, 0.24, 0.40];

    fig = figure( ...
        'Color', 'white', ...
        'Units', 'centimeters', ...
        'Position', [2, 2, 18, 9.5]);
    mainAxes = axes(fig, 'Position', [0.08, 0.20, 0.87, 0.74]);
    hold(mainAxes, 'on');

    drawStructuredBox(mainAxes, gridDimensions, physicalDimensions, ...
        paleBlue, gridBlue);

    wellTop = -0.35 * max(ly, lz);
    plot3(mainAxes, [injector(1), injector(1)], ...
        [injector(2), injector(2)], [injector(3), wellTop], '-', ...
        'Color', injectorBlue, 'LineWidth', 2.2);
    plot3(mainAxes, injector(1), injector(2), wellTop, '^', ...
        'MarkerSize', 7.5, 'MarkerFaceColor', injectorBlue, ...
        'MarkerEdgeColor', 'white', 'LineWidth', 0.8);

    plot3(mainAxes, [producer(1), producer(1)], ...
        [producer(2), producer(2)], [producer(3), wellTop], '-', ...
        'Color', producerBlue, 'LineWidth', 2.2);
    plot3(mainAxes, producer(1), producer(2), wellTop, 'v', ...
        'MarkerSize', 7.5, 'MarkerFaceColor', producerBlue, ...
        'MarkerEdgeColor', 'white', 'LineWidth', 0.8);

    xlabel(mainAxes, 'x (m)');
    ylabel(mainAxes, '');
    zlabel(mainAxes, 'z (m)');
    xlim(mainAxes, [-0.02 * lx, 1.02 * lx]);
    ylim(mainAxes, [-0.08 * ly, 1.08 * ly]);
    zlim(mainAxes, [wellTop - 0.08 * lz, 1.05 * lz]);
    xticks(mainAxes, linspace(0, lx, 7));
    yticks(mainAxes, linspace(0, ly, 3));
    zticks(mainAxes, linspace(0, lz, 3));
    view(mainAxes, -40, 25);
    camproj(mainAxes, 'orthographic');
    daspect(mainAxes, [1, 1, 1]);
    axis(mainAxes, 'vis3d');
    box(mainAxes, 'off');
    grid(mainAxes, 'off');
    set(mainAxes, ...
        'FontName', 'Arial', ...
        'FontSize', 8.5, ...
        'LineWidth', 0.7, ...
        'TickDir', 'out', ...
        'ZDir', 'reverse', ...
        'Color', 'white', ...
        'XColor', darkBlue, ...
        'YColor', darkBlue, ...
        'ZColor', darkBlue);
    annotation(fig, 'textbox', [0.012, 0.075, 0.10, 0.06], ...
        'String', 'y (m)', ...
        'LineStyle', 'none', ...
        'FontName', 'Arial', ...
        'FontSize', 8.5, ...
        'Color', darkBlue, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    camlight(mainAxes, 'headlight');

    stem = fullfile(outputDirectory, 'scw_kerogen_lmh_2d_geology');
    outputFiles = stem + [".png", ".pdf", ".svg", ".fig"];
    exportgraphics(fig, outputFiles(1), 'Resolution', 600, ...
        'BackgroundColor', 'white');
    exportgraphics(fig, outputFiles(2), 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    print(fig, outputFiles(3), '-dsvg', '-vector');
    savefig(fig, outputFiles(4));

    fprintf('Grid: %d x %d x %d cells, %.3f x %.3f x %.3f m\n', ...
        nx, ny, nz, lx, ly, lz);
    fprintf('Injector cell: i=%d, j=%d, k=0; producer cell: i=%d, j=%d, k=0\n', ...
        injectorI, primaryCentreJ, producerI, primaryCentreJ);
    fprintf('Wrote geological-model figures to %s\n', outputDirectory);
end

function value = readNumericConstant(sourceText, constantName)
    pattern = sprintf([ ...
        '(?m)^\\s*static\\s+constexpr\\s+(?:int|double)\\s+%s' ...
        '\\s*=\\s*([0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)\\s*;'], ...
        regexptranslate('escape', constantName));
    token = regexp(sourceText, pattern, 'tokens', 'once');
    if isempty(token)
        error('Could not read constant "%s" from the grid configuration.', ...
            constantName);
    end
    value = str2double(token{1});
end

function drawStructuredBox(ax, gridDimensions, physicalDimensions, ...
        faceColor, edgeColor)
    nx = gridDimensions(1);
    ny = gridDimensions(2);
    nz = gridDimensions(3);
    lx = physicalDimensions(1);
    ly = physicalDimensions(2);
    lz = physicalDimensions(3);

    xEdges = linspace(0, lx, nx + 1);
    yEdges = linspace(0, ly, ny + 1);
    zEdges = linspace(0, lz, nz + 1);

    [xTop, yTop] = meshgrid(xEdges, yEdges);
    addSurface(ax, xTop, yTop, zeros(size(xTop)), faceColor, edgeColor);
    addSurface(ax, xTop, yTop, lz * ones(size(xTop)), faceColor, edgeColor);

    [xSide, zSide] = meshgrid(xEdges, zEdges);
    addSurface(ax, xSide, zeros(size(xSide)), zSide, faceColor, edgeColor);
    addSurface(ax, xSide, ly * ones(size(xSide)), zSide, faceColor, edgeColor);

    [yEnd, zEnd] = meshgrid(yEdges, zEdges);
    addSurface(ax, zeros(size(yEnd)), yEnd, zEnd, faceColor, edgeColor);
    addSurface(ax, lx * ones(size(yEnd)), yEnd, zEnd, faceColor, edgeColor);
end

function surfaceHandle = addSurface(ax, x, y, z, faceColor, edgeColor)
    surfaceHandle = surf(ax, x, y, z, ...
        'FaceColor', faceColor, ...
        'EdgeColor', edgeColor, ...
        'FaceAlpha', 0.91, ...
        'LineWidth', 0.22, ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.72, ...
        'DiffuseStrength', 0.34, ...
        'SpecularStrength', 0.04);
end
