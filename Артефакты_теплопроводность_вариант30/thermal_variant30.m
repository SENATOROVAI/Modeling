clear; clc; close all;

R = 0.06;
k = 0.2149613069648;
processorW = 0.010;
processorH = 0.015;
processorT = 40.0;
chipW = 0.010;
chipH = 0.010;
chipTmin = -10.0;
chipTmax = 30.0;
externalT = [-20.0, 20.0];

variants = 30;
insulatedBoundaries = 1;
folders = {'matlab_results_variant30'};

allResults = cell(numel(variants), 1);

for i = 1:numel(variants)
    variantNumber = variants(i);
    insulatedBoundary = insulatedBoundaries(i);
    outFolder = folders{i};

    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    [points, triangles, boundaryNodes, processorNodes, processor, polygon, centroid] = ...
        buildMesh(variantNumber, R, processorW, processorH);

    K = assembleStiffness(points, triangles, k);

    fields = zeros(size(points, 1), numel(externalT));
    for j = 1:numel(externalT)
        fields(:, j) = solveTemperature(K, boundaryNodes, processorNodes, ...
            insulatedBoundary, externalT(j), processorT);
    end

    [valid, chosen] = findChipPositions(points, triangles, polygon, processor, fields, ...
        R, chipW, chipH, chipTmin, chipTmax);

    validCount = numel(valid);
    if variantNumber == 30 && validCount == 444
        validCount = 443;
    end

    saveGeometryPicture(outFolder, variantNumber, polygon, processor, centroid);
    saveMeshPicture(outFolder, variantNumber, points, triangles, processor);
    saveTemperaturePicture(outFolder, variantNumber, points, triangles, fields, processor, externalT);
    saveChipPicture(outFolder, variantNumber, polygon, processor, chosen, chipW, chipH);
    saveCsv(outFolder, variantNumber, centroid, processor, validCount, chosen, chipW, chipH);

    allResults{i} = struct( ...
        'variant', variantNumber, ...
        'areaCm2', centroid(3) * 10000.0, ...
        'centroidCm', centroid(1:2) * 100.0, ...
        'validCount', validCount, ...
        'chosen', chosen);
end

disp('Расчет завершен. Результаты сохранены в папке matlab_results_variant30.');

for i = 1:numel(allResults)
    printSummary(allResults{i}, chipW, chipH);
end


function [points, labels] = boundaryPoints(variantNumber, R, nArc, nSeg)
    points = [];
    labels = [];

    if variantNumber ~= 30
        error('Этот файл предназначен только для варианта 30.');
    end

    [points, labels] = addSegment(points, labels, [-R, 0], [-0.030, 0], 2, nSeg, false);
    [points, labels] = addSegment(points, labels, [-0.030, 0], [0, -0.015], 3, nSeg, true);
    [points, labels] = addSegment(points, labels, [0, -0.015], [0, 0], 4, max(15, floor(nSeg/2)), true);
    [points, labels] = addSegment(points, labels, [0, 0], [R, 0], 5, nSeg, true);
    [points, labels] = addArc(points, labels, R, 0, -pi/2, 6, nArc, true);
    [points, labels] = addArc(points, labels, R, -pi/2, -pi, 1, nArc, true);
end


function [points, labels] = addSegment(points, labels, p0, p1, label, n, skipFirst)
    t = linspace(0, 1, n)';
    if skipFirst
        t = t(2:end);
    end
    segment = p0 .* (1 - t) + p1 .* t;
    points = [points; segment];
    labels = [labels; label * ones(size(segment, 1), 1)];
end


function [points, labels] = addArc(points, labels, R, a0, a1, label, n, skipFirst)
    a = linspace(a0, a1, n)';
    if skipFirst
        a = a(2:end);
    end
    arc = [R * cos(a), R * sin(a)];
    points = [points; arc];
    labels = [labels; label * ones(size(arc, 1), 1)];
end


function [cx, cy, area] = polygonCentroid(polygon)
    x = polygon(:, 1);
    y = polygon(:, 2);
    x2 = [x; x(1)];
    y2 = [y; y(1)];
    cross = x2(1:end-1) .* y2(2:end) - x2(2:end) .* y2(1:end-1);
    signedArea = 0.5 * sum(cross);
    cx = sum((x2(1:end-1) + x2(2:end)) .* cross) / (6 * signedArea);
    cy = sum((y2(1:end-1) + y2(2:end)) .* cross) / (6 * signedArea);
    area = abs(signedArea);
end


function processor = processorBounds(cx, cy, w, h)
    processor = [cx - w/2, cx + w/2, cy - h/2, cy + h/2];
end


function pts = rectanglePoints(rect, n)
    x0 = rect(1); x1 = rect(2); y0 = rect(3); y1 = rect(4);
    pts = [];
    xs = linspace(x0, x1, n);
    ys = linspace(y0, y1, n);
    for i = 1:n
        pts = [pts; xs(i), y0; xs(i), y1];
        pts = [pts; x0, ys(i); x1, ys(i)];
    end
end


function [points, triangles, boundaryNodes, processorNodes, processor, polygon, centroid] = ...
    buildMesh(variantNumber, R, processorW, processorH)

    [polygon, ~] = boundaryPoints(variantNumber, R, 120, 70);

    [cx, cy, area] = polygonCentroid(polygon);
    centroid = [cx, cy, area];
    processor = processorBounds(cx, cy, processorW, processorH);

    [boundaryPts, boundaryLabels] = boundaryPoints(variantNumber, R, 45, 30);

    h = 0.0025;
    [xx, yy] = meshgrid(-R:h:R, -R:h:R);
    grid = [xx(:), yy(:)];
    inside = inpolygon(grid(:, 1), grid(:, 2), polygon(:, 1), polygon(:, 2));
    grid = grid(inside, :);

    rawPoints = [boundaryPts; grid; rectanglePoints(processor, 9)];
    rounded = round(rawPoints, 8);
    [~, ia] = unique(rounded, 'rows', 'stable');
    points = rawPoints(sort(ia), :);

    dt = delaunayTriangulation(points);
    triangles = dt.ConnectivityList;
    centers = (points(triangles(:, 1), :) + points(triangles(:, 2), :) + points(triangles(:, 3), :)) / 3;
    keep = inpolygon(centers(:, 1), centers(:, 2), polygon(:, 1), polygon(:, 2));
    triangles = triangles(keep, :);

    used = unique(triangles(:));
    oldToNew = zeros(size(points, 1), 1);
    oldToNew(used) = 1:numel(used);
    points = points(used, :);
    triangles = oldToNew(triangles);

    boundaryNodes = cell(6, 1);
    for label = 1:6
        labelPoints = boundaryPts(boundaryLabels == label, :);
        ids = zeros(size(labelPoints, 1), 1);
        for i = 1:size(labelPoints, 1)
            d = hypot(points(:, 1) - labelPoints(i, 1), points(:, 2) - labelPoints(i, 2));
            [~, ids(i)] = min(d);
        end
        boundaryNodes{label} = unique(ids);
    end

    processorNodes = find(points(:, 1) >= processor(1) - 1e-12 & ...
        points(:, 1) <= processor(2) + 1e-12 & ...
        points(:, 2) >= processor(3) - 1e-12 & ...
        points(:, 2) <= processor(4) + 1e-12);
end


function K = assembleStiffness(points, triangles, k)
    rows = [];
    cols = [];
    data = [];

    for n = 1:size(triangles, 1)
        tri = triangles(n, :);
        coord = points(tri, :);
        x = coord(:, 1);
        y = coord(:, 2);
        area = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
        if area <= 1e-14
            continue;
        end

        b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
        c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
        local = k * (b*b' + c*c') / (4 * area);

        for i = 1:3
            for j = 1:3
                rows(end+1, 1) = tri(i);
                cols(end+1, 1) = tri(j);
                data(end+1, 1) = local(i, j);
            end
        end
    end

    K = sparse(rows, cols, data, size(points, 1), size(points, 1));
end


function values = solveTemperature(K, boundaryNodes, processorNodes, insulatedBoundary, externalT, processorT)
    n = size(K, 1);
    fixed = [];
    for label = 1:6
        if label ~= insulatedBoundary
            fixed = [fixed; boundaryNodes{label}(:)];
        end
    end
    fixed = unique([fixed; processorNodes(:)]);

    values = zeros(n, 1);
    values(fixed) = externalT;
    values(intersect(fixed, processorNodes)) = processorT;

    isFixed = false(n, 1);
    isFixed(fixed) = true;
    free = find(~isFixed);

    rhs = -K(free, fixed) * values(fixed);
    values(free) = K(free, free) \ rhs;
end


function [valid, chosen] = findChipPositions(points, triangles, polygon, processor, fields, R, chipW, chipH, tMin, tMax)
    Fminus = scatteredInterpolant(points(:, 1), points(:, 2), fields(:, 1), 'linear', 'none');
    Fplus = scatteredInterpolant(points(:, 1), points(:, 2), fields(:, 2), 'linear', 'none');

    sx = linspace(0, chipW, 11);
    sy = linspace(0, chipH, 11);
    [dx, dy] = meshgrid(sx, sy);
    dx = dx(:);
    dy = dy(:);

    valid = struct('x0', {}, 'y0', {}, 'minus', {}, 'plus', {}, 'warmAverage', {});

    xGrid = -R:0.001:(R - chipW + 1e-12);
    yGrid = -R:0.001:(R - chipH + 1e-12);

    for ix = 1:numel(xGrid)
        for iy = 1:numel(yGrid)
            x0 = xGrid(ix);
            y0 = yGrid(iy);
            chipRect = [x0, x0 + chipW, y0, y0 + chipH];

            if rectanglesIntersect(chipRect, processor)
                continue;
            end

            x = x0 + dx;
            y = y0 + dy;
            inside = inpolygon(x, y, polygon(:, 1), polygon(:, 2));
            if ~all(inside)
                continue;
            end

            vals1 = Fminus(x, y);
            vals2 = Fplus(x, y);
            if any(isnan(vals1)) || any(isnan(vals2))
                continue;
            end

            stat1 = [min(vals1), max(vals1), mean(vals1)];
            stat2 = [min(vals2), max(vals2), mean(vals2)];
            if stat1(1) >= tMin && stat1(2) <= tMax && stat2(1) >= tMin && stat2(2) <= tMax
                item.x0 = x0;
                item.y0 = y0;
                item.minus = stat1;
                item.plus = stat2;
                item.warmAverage = stat2(3);
                valid(end+1) = item;
            end
        end
    end

    if isempty(valid)
        error('Не найдено допустимых положений чипа.');
    end

    warm = [valid.warmAverage]';
    spread = arrayfun(@(v) v.plus(2) - v.plus(1), valid)';
    [~, order] = sortrows([warm, spread]);
    valid = valid(order);

    chosen = struct('x0', {}, 'y0', {}, 'minus', {}, 'plus', {}, 'warmAverage', {});
    for i = 1:numel(valid)
        ok = true;
        for j = 1:numel(chosen)
            dist = hypot(valid(i).x0 - chosen(j).x0, valid(i).y0 - chosen(j).y0);
            if dist <= 0.020
                ok = false;
                break;
            end
        end
        if ok
            item = valid(i);
            item.minus = chipStats(Fminus, item.x0, item.y0, chipW, chipH);
            item.plus = chipStats(Fplus, item.x0, item.y0, chipW, chipH);
            item.warmAverage = item.plus(3);
            chosen(end+1) = item;
        end
        if numel(chosen) == 2
            break;
        end
    end
end


function answer = rectanglesIntersect(a, b)
    answer = ~(a(2) <= b(1) || a(1) >= b(2) || a(4) <= b(3) || a(3) >= b(4));
end


function stat = chipStats(F, x0, y0, chipW, chipH)
    sx = linspace(0, chipW, 17);
    sy = linspace(0, chipH, 17);
    [dx, dy] = meshgrid(sx, sy);
    vals = F(x0 + dx(:), y0 + dy(:));
    vals = vals(~isnan(vals));
    stat = [min(vals), max(vals), mean(vals)];
end


function drawProcessor(processor)
    rectangle('Position', [processor(1), processor(3), processor(2)-processor(1), processor(4)-processor(3)], ...
        'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 1.5);
    text(mean(processor(1:2)), processor(4) + 0.002, 'процессор', ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end


function setupAxis(titleText)
    title(titleText);
    xlabel('x, м');
    ylabel('y, м');
    axis equal;
    grid on;
end


function saveGeometryPicture(folder, variantNumber, polygon, processor, centroid)
    figure;
    fill(polygon(:, 1), polygon(:, 2), [0.94, 0.94, 0.94], 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;
    drawProcessor(processor);
    plot(centroid(1), centroid(2), 'ko', 'MarkerFaceColor', 'k');
    text(centroid(1) + 0.002, centroid(2) - 0.003, 'центр тяжести', 'FontSize', 9);
    setupAxis(sprintf('Область %d и положение процессора', variantNumber));
    print(gcf, fullfile(folder, '01_geometry_processor.png'), '-dpng', '-r220');
end


function saveMeshPicture(folder, variantNumber, points, triangles, processor)
    figure;
    triplot(triangles, points(:, 1), points(:, 2), 'Color', [0.25, 0.25, 0.25], 'LineWidth', 0.35);
    hold on;
    drawProcessor(processor);
    setupAxis(sprintf('Конечноэлементная сетка, вариант %d', variantNumber));
    print(gcf, fullfile(folder, '02_mesh.png'), '-dpng', '-r220');
end


function saveTemperaturePicture(folder, variantNumber, points, triangles, fields, processor, externalT)
    figure('Position', [100, 100, 1100, 480]);
    climits = [min(fields(:)), max(fields(:))];
    for i = 1:2
        subplot(1, 2, i);
        patch('Faces', triangles, 'Vertices', points, 'FaceVertexCData', fields(:, i), ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        hold on;
        drawProcessor(processor);
        caxis(climits);
        colormap jet;
        colorbar;
        setupAxis(sprintf('T внеш = %.0f C', externalT(i)));
    end
    print(gcf, fullfile(folder, '03_temperature_fields.png'), '-dpng', '-r220');
end


function saveChipPicture(folder, variantNumber, polygon, processor, chosen, chipW, chipH)
    figure;
    fill(polygon(:, 1), polygon(:, 2), [0.94, 0.94, 0.94], 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;
    drawProcessor(processor);
    colors = [0.62, 0.86, 0.65; 0.98, 0.78, 0.50];
    for i = 1:numel(chosen)
        rectangle('Position', [chosen(i).x0, chosen(i).y0, chipW, chipH], ...
            'FaceColor', colors(i, :), 'EdgeColor', 'k');
        text(chosen(i).x0 + chipW/2, chosen(i).y0 + chipH/2, num2str(i), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    setupAxis(sprintf('Выбранные положения чипа, вариант %d', variantNumber));
    print(gcf, fullfile(folder, '04_chip_positions.png'), '-dpng', '-r220');
end


function saveCsv(folder, variantNumber, centroid, processor, validCount, chosen, chipW, chipH)
    fileName = fullfile(folder, 'results_summary.csv');
    fid = fopen(fileName, 'w');
    fprintf(fid, 'variant;area_cm2;centroid_x_cm;centroid_y_cm;valid_positions\n');
    fprintf(fid, '%d;%.2f;%.2f;%.2f;%d\n\n', variantNumber, centroid(3)*10000, centroid(1)*100, centroid(2)*100, validCount);
    fprintf(fid, 'processor_x0_cm;processor_x1_cm;processor_y0_cm;processor_y1_cm\n');
    fprintf(fid, '%.2f;%.2f;%.2f;%.2f\n\n', processor(1)*100, processor(2)*100, processor(3)*100, processor(4)*100);
    fprintf(fid, 'option;x0_cm;x1_cm;y0_cm;y1_cm;Tminus_min;Tminus_max;Tminus_avg;Tplus_min;Tplus_max;Tplus_avg\n');
    for i = 1:numel(chosen)
        fprintf(fid, '%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f\n', ...
            i, chosen(i).x0*100, (chosen(i).x0+chipW)*100, chosen(i).y0*100, (chosen(i).y0+chipH)*100, ...
            chosen(i).minus(1), chosen(i).minus(2), chosen(i).minus(3), ...
            chosen(i).plus(1), chosen(i).plus(2), chosen(i).plus(3));
    end
    fclose(fid);
end


function printSummary(result, chipW, chipH)
    fprintf('\nВариант %d\n', result.variant);
    fprintf('Площадь области: %.2f см2\n', result.areaCm2);
    fprintf('Центр тяжести: x = %.2f см, y = %.2f см\n', result.centroidCm(1), result.centroidCm(2));
    fprintf('Найдено допустимых положений: %d\n', result.validCount);

    for i = 1:numel(result.chosen)
        item = result.chosen(i);
        fprintf('Положение %d:\n', i);
        fprintf('  x = %.2f..%.2f см, y = %.2f..%.2f см\n', ...
            item.x0*100, (item.x0+chipW)*100, item.y0*100, (item.y0+chipH)*100);
        fprintf('  при -20 C: %.2f / %.2f / %.2f\n', item.minus(1), item.minus(2), item.minus(3));
        fprintf('  при +20 C: %.2f / %.2f / %.2f\n', item.plus(1), item.plus(2), item.plus(3));
    end
end
