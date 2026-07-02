function scene = planner_fangzhen_scene(repoRoot)
% 供 PRMO_fangzhen / planner_travel_*_fangzhen 使用：OSM 建筑物 + 禁飞盒 + 人流圆柱 + order 环境用于目标3
    persistent cachedRoot cachedScene
    if nargin < 1 || isempty(repoRoot)
        repoRoot = getenv('PLATEMO_FANGZHEN_ROOT');
    end
    if isempty(repoRoot)
        repoRoot = 'D:\PlatEMO-master-using';
    end
    if ~isempty(cachedScene) && ~isempty(cachedRoot) && strcmp(cachedRoot, repoRoot)
        scene = cachedScene;
        return;
    end
    orderFile = fullfile(repoRoot, 'PlatEMO', 'forOrderNew26', 'Order_Map', 'order_data_fangzhen.m');
    if ~isfile(orderFile)
        error('找不到 %s', orderFile);
    end
    eval(fileread(orderFile)); %#ok<EVLDIR>
    bmat = getenv('PLATEMO_FANGZHEN_BMAT');
    if isempty(bmat)
        bmat = fullfile(repoRoot, 'buildings_from_osm.mat');
        if ~isfile(bmat)
            pl = fullfile(repoRoot, 'PlatEMO');
            alt = fullfile(pl, 'buildings_from_osm.mat');
            if isfile(alt)
                bmat = alt;
            end
        end
    end
    if ~isfile(bmat)
        alt2 = fullfile(repoRoot, '26article', 'map', 'buildings_from_osm.mat');
        if isfile(alt2)
            bmat = alt2;
        else
            error('找不到建筑物 mat。请将 buildings_from_osm.mat 放在工程根、PlatEMO 下，或 setenv(''PLATEMO_FANGZHEN_BMAT'', 完整路径)');
        end
    end
    Sb = load(bmat, 'buildings');
    buildings = Sb.buildings;
    % 禁飞：长方体（硬避障，追加到 buildings）
    for iz = 1:size(forbidden_zones, 1)
        r = forbidden_zones(iz, :);
        buildings = [buildings; fangzhen_rect8_row(r(1), r(2), r(3), r(4), 0, 2400)]; %#ok<AGROW>
    end
    scene.buildings = buildings;
    scene.numBuildings = size(buildings, 1);
    % 人流区：按圆柱筛（与 checknode 中圆柱一致，高度取 z 上限）
    nc = size(crowded_zones, 1);
    scene.cylinders = zeros(max(0, nc), 4);
    for ic = 1:nc
        scene.cylinders(ic, :) = [crowded_zones(ic, 1:2), crowded_zones(ic, 3), 2400];
    end
    scene.numCylinders = size(scene.cylinders, 1);
    scene.pyramids = zeros(0, 24);
    scene.numPyramids = 0;
    scene.spheres = zeros(0, 4);
    scene.numSpheres = 0;
    % jinfei.mat 在 checknode 中只读第 1 个长方体，用第一块禁飞区
    r0 = forbidden_zones(1, :);
    scene.jinfei = fangzhen_rect8_row(r0(1), r0(2), r0(3), r0(4), 0, 2400);
    scene.numjinfei = 1;
    Rrect = [forbidden_zones; noise_zones];
    scene.xLimits = [Rrect(:, 1), Rrect(:, 3)];
    scene.yLimits = [Rrect(:, 2), Rrect(:, 4)];
    scene.circleCenter = crowded_zones(:, 1:2);
    scene.radius = crowded_zones(:, 3);

    cachedRoot = repoRoot;
    cachedScene = scene;
end
