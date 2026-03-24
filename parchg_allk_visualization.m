% =========================================================================
% VASP PARCHG / ALLK 3D Visualization Script
% 从结构头读取晶格与原子，从 ALLK 数据块读取三维电荷密度矩阵
% 改进点：
%   1) 原子显示与电荷密度裁剪彻底解耦，避免 Ti/O 因裁剪消失
%   2) X 方向采用“宽窗口+边缘剔除”策略：保留中间 Ti-Ti 成键电荷，删除最左/最右边缘噪声
% =========================================================================
clear; clc; close all;

%% 1. 核心可调参数区 (User Configurations)
filename = 'PARCHG.0065.ALLK';

% --- 图像与等值面设置 ---
isovalue = 5.0;
interp_factor = 2;
iso_color = [1.0 1.0 0.1];
iso_alpha = 0.6;

% --- 原子显示范围（尽量完整显示）---
% 使用 2x 视场，确保 Ti/O 不因边界条件漏掉
atom_x = [-0.5, 1.5];
atom_y = [-0.1, 1.1];
atom_z = [-0.1, 1.1];

% --- 等值面采样范围（先宽后裁）---
% 先覆盖与原子相同的大窗口，再通过 iso_trim_x 删除左右边缘
iso_x = atom_x;
iso_y = [0, 1];
iso_z = [0, 1];

% 仅删除 X 方向最左/最右边缘（保留中间 Ti-Ti 连续电荷）
enable_iso_edge_trim = true;
iso_trim_x = [-0.35, 1.35];

% --- 空间掩膜（可选）---
enable_spatial_mask = false;
mask_x = [0.0, 1.0];
mask_y = [0.2, 0.8];
mask_z = [0.0, 1.0];

% --- 原子显示设置 ---
show_atoms = true;
atom_colors = [
    0.6 0.6 0.8;  % Ti
    0.9 0.2 0.2   % O
];
atom_radii = [0.6, 0.4];

% --- 几何成键(Bonds)设置 ---
show_bonds = true;
bond_cutoff = 2.8;
bond_radius = 0.08;

%% 2. 解析数据文件
fprintf('正在读取文件: %s\n', filename);
fid = fopen(filename, 'r');
assert(fid > 0, '无法打开文件: %s', filename);

sys_name = strtrim(fgetl(fid)); %#ok<NASGU>
scale = str2double(fgetl(fid));

lat = zeros(3,3);
for i = 1:3
    lat(i,:) = str2num(fgetl(fid)) * scale; %#ok<ST2NM>
end

elem_line = strtrim(fgetl(fid));
if isstrprop(elem_line(1), 'alpha')
    counts = str2num(fgetl(fid)); %#ok<ST2NM>
else
    counts = str2num(elem_line); %#ok<ST2NM>
end
nAtoms = sum(counts);

atom_types = zeros(nAtoms, 1);
idx = 1;
for i = 1:length(counts)
    atom_types(idx : idx+counts(i)-1) = i;
    idx = idx + counts(i);
end

mode_line = strtrim(fgetl(fid));
if strcmpi(mode_line(1), 'S')
    mode_line = strtrim(fgetl(fid));
end
is_direct = strcmpi(mode_line(1), 'D');

coords = zeros(nAtoms, 3);
for i = 1:nAtoms
    temp = str2num(fgetl(fid)); %#ok<ST2NM>
    coords(i,:) = temp(1:3);
end
if ~is_direct
    coords = coords / lat;
end

next_line = strtrim(fgetl(fid));
while isempty(next_line)
    next_line = strtrim(fgetl(fid));
end
grid_dim = str2num(next_line); %#ok<ST2NM>
nx = grid_dim(1); ny = grid_dim(2); nz = grid_dim(3);

data_raw = fscanf(fid, '%f');
fclose(fid);

chg = reshape(data_raw(1:nx*ny*nz), [nx, ny, nz]);
fprintf('解析完毕。共 %d 个原子，网格: %d x %d x %d\n', nAtoms, nx, ny, nz);

%% 3. 周期扩充与插值（X 方向扩到 3 倍以支持 [-0.5,1.5] 视场）
fprintf('正在进行数据网格插值...\n');

chg_3x = repmat(chg, [3, 1, 1]);
x_base = (0:nx-1)' / nx;
x_3x = [x_base - 1; x_base; x_base + 1];

chg_pad = padarray(chg_3x, [0, 2, 2], 'circular', 'both');
y_base = (0:ny-1)' / ny;
y_pad = [y_base(end-1:end)-1; y_base; y_base(1:2)+1];
z_base = (0:nz-1)' / nz;
z_pad = [z_base(end-1:end)-1; z_base; z_base(1:2)+1];

margin = 0.1;
x_idx = find(x_3x >= iso_x(1)-margin & x_3x <= iso_x(2)+margin);
x_target = x_3x(x_idx);
chg_target = chg_pad(x_idx, :, :);

[Xg, Yg, Zg] = ndgrid(x_target, y_pad, z_pad);

n_xq = max(2, round((iso_x(2)-iso_x(1))*nx) * interp_factor);
n_yq = max(2, round((iso_y(2)-iso_y(1))*ny) * interp_factor);
n_zq = max(2, round((iso_z(2)-iso_z(1))*nz) * interp_factor);

xq = linspace(iso_x(1), iso_x(2), n_xq);
yq = linspace(iso_y(1), iso_y(2), n_yq);
zq = linspace(iso_z(1), iso_z(2), n_zq);
[Xq, Yq, Zq] = ndgrid(xq, yq, zq);

chg_fine = interpn(Xg, Yg, Zg, chg_target, Xq, Yq, Zq, 'spline');

%% 4. 空间掩膜 / 裁剪
if enable_spatial_mask
    mask = (Xq >= mask_x(1)) & (Xq <= mask_x(2)) & ...
           (Yq >= mask_y(1)) & (Yq <= mask_y(2)) & ...
           (Zq >= mask_z(1)) & (Zq <= mask_z(2));
    chg_fine(~mask) = -1e9;
end

if enable_iso_edge_trim
    edge_keep = (Xq >= iso_trim_x(1)) & (Xq <= iso_trim_x(2));
    chg_fine(~edge_keep) = -1e9;
end

%% 5. 渲染
fprintf('正在渲染三维图像与原子...\n');
figure('Name', 'PARCHG/ALLK Visualization', 'Color', 'w', 'Position', [100, 100, 900, 700]);
hold on;

fv = isosurface(Xq, Yq, Zq, chg_fine, isovalue);
if ~isempty(fv.vertices)
    fv.vertices = fv.vertices * lat;
    p = patch(fv);
    p.FaceColor = iso_color;
    p.EdgeColor = 'none';
    p.FaceAlpha = iso_alpha;
    p.AmbientStrength = 0.6;
    p.DiffuseStrength = 0.8;
end

plot_cell_box(lat, [0, 1], [0, 1], [0, 1], '-', 'k');

%% 6. 收集并绘制原子与成键
if show_atoms
    % 扩展平移范围，确保边界 Ti/O 不漏
    [dx, dy, dz] = ndgrid(-1:2, -1:2, -1:2);
    shifts = [dx(:), dy(:), dz(:)];

    vis_coords = zeros(0,3);
    vis_types = zeros(0,1);
    eps_b = 0.03;

    for i = 1:nAtoms
        for j = 1:size(shifts, 1)
            new_frac = coords(i,:) + shifts(j,:);
            if new_frac(1) >= atom_x(1)-eps_b && new_frac(1) <= atom_x(2)+eps_b && ...
               new_frac(2) >= atom_y(1)-eps_b && new_frac(2) <= atom_y(2)+eps_b && ...
               new_frac(3) >= atom_z(1)-eps_b && new_frac(3) <= atom_z(2)+eps_b
                vis_coords(end+1,:) = new_frac * lat; %#ok<AGROW>
                vis_types(end+1,1) = atom_types(i); %#ok<AGROW>
            end
        end
    end

    if show_bonds
        nVis = size(vis_coords, 1);
        for k1 = 1:nVis
            for k2 = k1+1:nVis
                dist = norm(vis_coords(k1,:) - vis_coords(k2,:));
                if dist > 0.1 && dist <= bond_cutoff
                    c1 = atom_colors(vis_types(k1), :);
                    c2 = atom_colors(vis_types(k2), :);
                    draw_bicolor_bond(vis_coords(k1,:), vis_coords(k2,:), bond_radius, c1, c2);
                end
            end
        end
    end

    [sx, sy, sz] = sphere(30);
    for k = 1:size(vis_coords, 1)
        r = atom_radii(vis_types(k));
        c = atom_colors(vis_types(k), :);
        cart = vis_coords(k, :);
        surf(sx*r + cart(1), sy*r + cart(2), sz*r + cart(3), ...
             'FaceColor', c, 'EdgeColor', 'none', ...
             'AmbientStrength', 0.6, 'DiffuseStrength', 0.8, ...
             'SpecularStrength', 0.3);
    end
end

%% 7. 光照与视角
view(3);
axis equal;
axis off;
camlight('headlight');
camlight('right');
lighting gouraud;
material default;
title(sprintf('Isovalue: %g', isovalue), 'Interpreter', 'none', 'FontSize', 12);
fprintf('执行完毕。\n');

% =========================================================================
function plot_cell_box(lat, xlims, ylims, zlims, line_style, line_color)
    x = [xlims(1) xlims(2)];
    y = [ylims(1) ylims(2)];
    z = [zlims(1) zlims(2)];
    [X, Y, Z] = ndgrid(x, y, z);
    frac_corners = [X(:), Y(:), Z(:)];
    cart_corners = frac_corners * lat;

    edges = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; 4 8; 5 6; 5 7; 6 8; 7 8];
    for i = 1:size(edges, 1)
        p1 = cart_corners(edges(i,1), :);
        p2 = cart_corners(edges(i,2), :);
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
            'Color', line_color, 'LineStyle', line_style, 'LineWidth', 1.0);
    end
end

% =========================================================================
function draw_bicolor_bond(p1, p2, r, color1, color2)
    v = p2 - p1;
    L = norm(v);
    dir = v / L;

    [xc, yc, zc] = cylinder(r, 20);

    z_axis = [0, 0, 1];
    axis_rot = cross(z_axis, dir);
    angle = acos(max(-1,min(1,dot(z_axis, dir))));

    if norm(axis_rot) < 1e-6
        if dot(z_axis, dir) < 0
            R = [1 0 0; 0 -1 0; 0 0 -1];
        else
            R = eye(3);
        end
    else
        axis_rot = axis_rot / norm(axis_rot);
        K = [0 -axis_rot(3) axis_rot(2); axis_rot(3) 0 -axis_rot(1); -axis_rot(2) axis_rot(1) 0];
        R = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
    end

    zc1 = zc * (L/2);
    pts1 = [xc(:) yc(:) zc1(:)] * R';
    X1 = reshape(pts1(:,1), size(xc)) + p1(1);
    Y1 = reshape(pts1(:,2), size(yc)) + p1(2);
    Z1 = reshape(pts1(:,3), size(zc)) + p1(3);
    surf(X1, Y1, Z1, 'FaceColor', color1, 'EdgeColor', 'none', 'AmbientStrength', 0.6, 'DiffuseStrength', 0.8);

    zc2 = zc * (L/2) + (L/2);
    pts2 = [xc(:) yc(:) zc2(:)] * R';
    X2 = reshape(pts2(:,1), size(xc)) + p1(1);
    Y2 = reshape(pts2(:,2), size(yc)) + p1(2);
    Z2 = reshape(pts2(:,3), size(zc)) + p1(3);
    surf(X2, Y2, Z2, 'FaceColor', color2, 'EdgeColor', 'none', 'AmbientStrength', 0.6, 'DiffuseStrength', 0.8);
end
