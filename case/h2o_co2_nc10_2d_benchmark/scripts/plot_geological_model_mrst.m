function plot_geological_model_mrst(outputDirectory, mrstRoot)
%PLOT_GEOLOGICAL_MODEL_MRST 绘制 H2O-CO2-nC10 二维算例地质模型。
%   使用 MRST 的 cartGrid、computeGeometry 与 plotGrid 重建正式算例的
%   60x20x1 网格，并输出适合论文排版的 PNG、PDF、SVG 和 MATLAB FIG。

    arguments
        outputDirectory (1, 1) string = ""
        mrstRoot (1, 1) string = ""
    end

    scriptDirectory = fileparts(mfilename('fullpath'));
    caseDirectory = fileparts(scriptDirectory);
    repositoryRoot = fileparts(fileparts(caseDirectory));

    if strlength(outputDirectory) == 0
        outputDirectory = fullfile(repositoryRoot, 'outputs', 'figures', 'geological_model');
    end
    if strlength(mrstRoot) == 0
        mrstRoot = string(getenv('MRST_ROOT'));
    end
    if strlength(mrstRoot) == 0
        mrstRoot = "D:\DOCUMENTS\研究生\MRST\mrst-2024b";
    end
    if ~isfolder(mrstRoot)
        error('MRST root does not exist: %s', mrstRoot);
    end
    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    addpath(mrstRoot, '-begin');
    if exist('cartGrid', 'file') ~= 2
        feval('startup');
    end

    gridDimensions = [60, 20, 1];
    physicalDimensions = [300, 100, 5];
    G = cartGrid(gridDimensions, physicalDimensions);
    G = computeGeometry(G);

    injectorCell = sub2ind(G.cartDims, 1, 10, 1);
    producerCell = sub2ind(G.cartDims, 60, 10, 1);
    injector = G.cells.centroids(injectorCell, :);
    producer = G.cells.centroids(producerCell, :);

    paleBlue = [0.78, 0.87, 0.94];
    gridBlue = [0.43, 0.59, 0.72];
    darkBlue = [0.12, 0.31, 0.48];
    injectorBlue = [0.16, 0.55, 0.72];
    producerBlue = [0.08, 0.24, 0.40];

    fig = figure( ...
        'Color', 'white', ...
        'Units', 'centimeters', ...
        'Position', [2, 2, 17, 9.5], ...
        'Renderer', 'painters');

    mainAxes = axes(fig, 'Position', [0.07, 0.13, 0.89, 0.82]);
    hold(mainAxes, 'on');
    axes(mainAxes); %#ok<LAXES>
    gridHandle = plotGrid( ...
        G, ...
        'FaceColor', paleBlue, ...
        'EdgeColor', gridBlue, ...
        'FaceAlpha', 0.91, ...
        'LineWidth', 0.22);
    if isgraphics(gridHandle)
        set(gridHandle, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.72, ...
            'DiffuseStrength', 0.34, 'SpecularStrength', 0.04);
    end

    wellTop = -10.0;
    plot3(mainAxes, [injector(1), injector(1)], [injector(2), injector(2)], ...
        [injector(3), wellTop], '-', 'Color', injectorBlue, 'LineWidth', 2.2);
    plot3(mainAxes, injector(1), injector(2), wellTop, '^', ...
        'MarkerSize', 7.5, 'MarkerFaceColor', injectorBlue, ...
        'MarkerEdgeColor', 'white', 'LineWidth', 0.8);
    plot3(mainAxes, [producer(1), producer(1)], [producer(2), producer(2)], ...
        [producer(3), wellTop], '-', 'Color', producerBlue, 'LineWidth', 2.2);
    plot3(mainAxes, producer(1), producer(2), wellTop, 'v', ...
        'MarkerSize', 7.5, 'MarkerFaceColor', producerBlue, ...
        'MarkerEdgeColor', 'white', 'LineWidth', 0.8);

    xlabel(mainAxes, 'x (m)');
    ylabel(mainAxes, 'y (m)');
    zlabel(mainAxes, 'z (m)');
    xlim(mainAxes, [-5, 305]);
    ylim(mainAxes, [-5, 105]);
    zlim(mainAxes, [-13, 6]);
    xticks(mainAxes, [0, 100, 200, 300]);
    yticks(mainAxes, [0, 50, 100]);
    zticks(mainAxes, [0, 5]);
    view(mainAxes, -40, 27);
    camproj(mainAxes, 'orthographic');
    daspect(mainAxes, [1, 1, 0.25]);
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
    camlight(mainAxes, 'headlight');

    stem = fullfile(outputDirectory, '28_geological_model_mrst');
    exportgraphics(fig, stem + ".png", 'Resolution', 600, ...
        'BackgroundColor', 'white');
    exportgraphics(fig, stem + ".pdf", 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    print(fig, stem + ".svg", '-dsvg', '-painters');
    savefig(fig, stem + ".fig");
    close(fig);

    fprintf('Wrote MRST geological-model figure to %s\n', outputDirectory);
end
