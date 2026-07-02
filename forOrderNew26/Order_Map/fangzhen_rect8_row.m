function row = fangzhen_rect8_row(xmin, ymin, xmax, ymax, z_bot, z_top)
% 与 import_bbbike_buildings_to_mat 一致的底面/顶面 8 顶点 1x24
v1 = [xmin; ymin; z_bot];
v2 = [xmax; ymin; z_bot];
v3 = [xmax; ymax; z_bot];
v4 = [xmin; ymax; z_bot];
v5 = [xmin; ymin; z_top];
v6 = [xmax; ymin; z_top];
v7 = [xmax; ymax; z_top];
v8 = [xmin; ymax; z_top];
V = [v1, v2, v3, v4, v5, v6, v7, v8];
row = reshape(V, 1, 24);
end
