function p = fangzhen_id_to_dot(i0, dotss, zf)
%FANGZHEN_ID_TO_DOT  节点 id -> [x y z]（与 order_data_fangzhen 的 dotss 一致）
% i0: 0 = 第 11 行(配送中心), 1~10 = 第 1~10 行
    if i0 < 0 || i0 > 10
        error('节点编号必须在 0~10, 得 %d', i0);
    end
    if i0 == 0
        row = 11;
    else
        row = i0;
    end
    p = [dotss(row, 1:2), zf];
end
