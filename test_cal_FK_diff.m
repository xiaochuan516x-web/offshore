function WaveExcitationForce = test_cal_FK_diff(Global, Time, Farm, SampleTime, Wave, WaveDataMode, vel_history)
% test_cal_FK_diff - 水动力计算
%
% 计算逻辑：
% 1) STL 必须包含全部部件；ActiveComponents 是唯一选择器。
% 2) FK：逐部件计算，只累加 ActiveComponents。
% 3) Viscous drag：逐部件计算，只累加 ActiveComponents。
% 4) Radiation：AQWA/K 对完整平台整体计算，再按"全部部件"的瞬时 WetArea 分配，
%    最后只累加 ActiveComponents。
% 5) Diffraction：AQWA 对完整平台整体计算，再按"全部部件"的瞬时 WetArea 分配，
%    最后只累加 ActiveComponents。
% 6) 不再使用 Enabled / EnableFK / EnableDrag / EnableRadiation / EnableDiffraction。

%% 1. 固定风机数量和输出
num_turbines = 1;
WaveExcitationForce = zeros(2,num_turbines);

%% 2. Persistent 数据
persistent Turb K_matrix diffData last_diag_time ...
    all_Area all_Center all_Tnorm all_V1 all_V2 all_V3 all_Volume ...
    comp_start_idx comp_end_idx num_components ActiveComponents

if isempty(Turb)
    Pars = load('Parameters_Simu.mat');
    Turb = Pars.Turb;
    K_matrix = Pars.Turb(1).K;
    diffData = Pars.Turb(1).diffData;
    ActiveComponents = Pars.ActiveComponents;

    Mesh = load('STLMeshData.mat');
    all_Area = Mesh.all_Area;
    all_Center = Mesh.all_Center;
    all_Tnorm = Mesh.all_Tnorm;
    all_Volume = Mesh.all_Volume;
    comp_start_idx = Mesh.comp_start_idx;
    comp_end_idx = Mesh.comp_end_idx;
    num_components = Mesh.num_parts;

    all_V1 = Mesh.all_V1;
    all_V2 = Mesh.all_V2;
    all_V3 = Mesh.all_V3;
end

%% 3. 整理 ActiveComponents
ActiveComponents = ActiveComponents(:)';
active_list = zeros(1,num_components);
nActive = 0;
for i = 1:length(ActiveComponents)
    comp_id = ActiveComponents(i);
    if comp_id >= 1 && comp_id <= num_components && comp_id == floor(comp_id)
        duplicate = false;
        for j = 1:nActive
            if active_list(j) == comp_id
                duplicate = true;
                break;
            end
        end
        if ~duplicate
            nActive = nActive + 1;
            active_list(nActive) = comp_id;
        end
    end
end
if nActive == 0
    return;
end

%% 4. Wave 数据
kx = Wave.kx(:)';
ky = Wave.ky(:)';
omega = Wave.omega(:)';

if WaveDataMode == 0
    amp = Wave.amp0(:)';
    phase = Wave.phase0(:)' - Wave.omega(:)' * Time;
else
    if Time == 0
        amp = Wave.amp0(:)';
        phase = Wave.phase0(:)';
    else
        current_step = int32(floor(Time/SampleTime + 1e-5));
        filename = sprintf('%swave_%04d.dat',Wave.wave_case,current_step);
        fileID = fopen(filename,'r');
        if fileID == -1
            amp = Wave.amp0(:)';
            phase = Wave.phase0(:)';
        else
            data = fscanf(fileID,'%f');
            fclose(fileID);
            if isempty(data)
                amp = Wave.amp0(:)';
                phase = Wave.phase0(:)';
            else
                temp_matrix = reshape(data,2,[]);
                amp = temp_matrix(1,:);
                phase = temp_matrix(2,:);
            end
        end
    end
end

nWave = min([length(kx),length(ky),length(amp),length(phase),length(omega)]);
if nWave <= 0
    return;
end
kx = kx(1:nWave);
ky = ky(1:nWave);
omega = omega(1:nWave);
amp = amp(1:nWave);
phase = phase(1:nWave);
wave_dat = [kx;ky;amp;phase];

%% 5. 水动力
for TurbNum = 1:num_turbines
    dof = zeros(2,1);
    dof(1) = Turb(TurbNum).PosVec(1);
    dof(2) = Turb(TurbNum).PosVec(2);

    rho = Global.rho;
    g = Global.g;
    cg = Turb(TurbNum).NeutralPosVec;

    % 当前瞬时波面：仍按原来的 center + 当前平面位置计算。
    elv = waveElevation(dof,all_Center,cg,amp,phase,kx,ky,0);

    % 关键：WetArea 对"全部部件"计算，不能因为 ActiveComponents 而删除未选部件。
    if ~isempty(all_V1)
        WetArea = calculate_dynamic_wet_area_exact(...
            all_Area,all_V1,all_V2,all_V3,dof,cg,amp,phase,kx,ky,...
            comp_start_idx,comp_end_idx,num_components);
    else
        WetArea = calculate_dynamic_wet_area_center(...
            all_Area,all_Center,elv,comp_start_idx,comp_end_idx,num_components);
    end

    FK_force = calculate_FK_force_by_component(...
        active_list,nActive,all_Area,all_Center,all_Tnorm,...
        comp_start_idx,comp_end_idx,num_components,elv,rho,g,amp,phase,kx,ky,dof,cg);

    Viscous_force = calculate_viscous_force_by_component(...
        active_list,nActive,TurbNum,all_Center,comp_start_idx,comp_end_idx,...
        num_components,WetArea,Turb(TurbNum).VelVec,rho,g,amp,phase,kx,ky,omega,Turb,dof,cg);

    % Radiation/Diffraction 先得到完整平台整体结果，再按全部 WetArea 分配。
    Rad_force = calculate_radiation_force_by_wet_area(...
        WetArea,active_list,nActive,vel_history,K_matrix,SampleTime);

    Diff_force = calculate_diffraction_force_by_wet_area(...
        wave_dat,diffData,dof,g,omega,WetArea,active_list,nActive);

    WaveExcitationForce(:,TurbNum) = FK_force + Viscous_force + Rad_force + Diff_force;

    % Continuous-time diagnostic. For a 20 s regular wave, sampling only at
    % 10 s intervals lands at half-period increments and naturally shows
    % alternating signs. This logger samples every 1 s and separates terms.
    if isempty(last_diag_time)
        last_diag_time = -inf;
    end
    if Time >= 0 && (Time - last_diag_time >= 1.0 - 1e-9)
        fprintf(['[HydroDiag] T=%8.3f | FK=[% .4e,% .4e] | Visc=[% .4e,% .4e] ' ...
                 '| Rad=[% .4e,% .4e] | Diff=[% .4e,% .4e] | Total=[% .4e,% .4e]\n'], ...
            Time, FK_force(1), FK_force(2), Viscous_force(1), Viscous_force(2), ...
            Rad_force(1), Rad_force(2), Diff_force(1), Diff_force(2), ...
            WaveExcitationForce(1,TurbNum), WaveExcitationForce(2,TurbNum));
        last_diag_time = Time;
    end
end
end

%% ========================================================================
% 子函数1：波面高度
% ========================================================================
function f = waveElevation(x,center,cg,AH1,kphase,kx,ky,t)
    [f,~] = calc_elev_wang(x,center,cg,AH1,kphase,kx,ky,t);
end

function [f,sf] = calc_elev_wang(x,center,cg,AH1,kphase,kx,ky,t)
    center = offsetXYZ(center,x);
    center = offsetXYZ(center,cg);
    [f,sf] = waveElev_wang(center,AH1,kx,ky,kphase,t);
end

function [f,sf] = waveElev_wang(center,AH1,kx,ky,kphase,t)
    f = zeros(size(center,1),1);
    sf = zeros(size(center,1),length(kx));
    for ic = 1:size(center,1)
        for ik = 1:length(kx)
            sf(ic,ik) = AH1(ik)*cos(kx(ik)*center(ic,1) + ...
                ky(ik)*center(ic,2) + kphase(ik));
        end
    end
    f = sum(sf,2);
end

function verts = offsetXYZ(verts,x)
    verts(:,1) = verts(:,1) + x(1);
    verts(:,2) = verts(:,2) + x(2);
end

%% ========================================================================
% 子函数2：动态 WetArea
% 三角形水下面积：把每个三角形按三个顶点的 signed depth 做线性裁剪。
% 水面在三角形内部采用三个顶点波面值的线性插值。
% ========================================================================
function WetArea = calculate_dynamic_wet_area_exact(all_Area,all_V1,all_V2,all_V3,...
    dof,cg,amp,phase,kx,ky,comp_start_idx,comp_end_idx,num_components)

    WetArea = zeros(num_components,1);
    nTri = length(all_Area);

    for i = 1:nTri
        p1 = all_V1(i,:) + [dof(1)+cg(1),dof(2)+cg(2),0];
        p2 = all_V2(i,:) + [dof(1)+cg(1),dof(2)+cg(2),0];
        p3 = all_V3(i,:) + [dof(1)+cg(1),dof(2)+cg(2),0];

        e1 = p1(3) - local_wave_elevation_xy(p1(1),p1(2),amp,phase,kx,ky);
        e2 = p2(3) - local_wave_elevation_xy(p2(1),p2(2),amp,phase,kx,ky);
        e3 = p3(3) - local_wave_elevation_xy(p3(1),p3(2),amp,phase,kx,ky);

        if e1 <= 0 && e2 <= 0 && e3 <= 0
            wet = all_Area(i);
        elseif e1 > 0 && e2 > 0 && e3 > 0
            wet = 0;
        else
            wet = clipped_triangle_area(p1,p2,p3,e1,e2,e3);
        end

        comp = find_component_from_index(i,comp_start_idx,comp_end_idx,num_components);
        if comp >= 1
            WetArea(comp) = WetArea(comp) + wet;
        end
    end
end

function WetArea = calculate_dynamic_wet_area_center(all_Area,all_Center,elv,...
    comp_start_idx,comp_end_idx,num_components)
    WetArea = zeros(num_components,1);
    for comp = 1:num_components
        s = comp_start_idx(comp);
        e = comp_end_idx(comp);
        if s <= e
            for i = s:e
                if all_Center(i,3) <= elv(i)
                    WetArea(comp) = WetArea(comp) + all_Area(i);
                end
            end
        end
    end
end

function zeta = local_wave_elevation_xy(x,y,amp,phase,kx,ky)
    zeta = 0;
    for ik = 1:length(kx)
        zeta = zeta + amp(ik)*cos(kx(ik)*x + ky(ik)*y + phase(ik));
    end
end

function wet = clipped_triangle_area(p1,p2,p3,d1,d2,d3)
    % 把三角形裁剪成 signed depth <= 0 的多边形，再用叉积求面积。
    pts = zeros(6,3);
    n = 0;
    P = [p1;p2;p3];
    D = [d1;d2;d3];
    for i = 1:3
        j = i + 1;
        if j > 3, j = 1; end
        Pi = P(i,:); Pj = P(j,:);
        Di = D(i); Dj = D(j);
        if Di <= 0
            n = n + 1;
            pts(n,:) = Pi;
        end
        if (Di < 0 && Dj > 0) || (Di > 0 && Dj < 0)
            alpha = Di/(Di-Dj);
            n = n + 1;
            pts(n,:) = Pi + alpha*(Pj-Pi);
        elseif Di == 0 && Dj > 0
            % Di=0 的点已经加入，不重复添加。
        end
    end
    if n < 3
        wet = 0;
        return;
    end
    Q = pts(1:n,:);
    wet = polygon_area_3d(Q);
end

function A = polygon_area_3d(P)
    A = 0;
    if size(P,1) < 3, return; end
    p0 = P(1,:);
    for i = 2:size(P,1)-1
        A = A + 0.5*norm(cross(P(i,:)-p0,P(i+1,:)-p0));
    end
end

function comp = find_component_from_index(idx,comp_start_idx,comp_end_idx,num_components)
    comp = 0;
    for k = 1:num_components
        if idx >= comp_start_idx(k) && idx <= comp_end_idx(k)
            comp = k;
            return;
        end
    end
end

%% ========================================================================
% 子函数3：FK 力
% ========================================================================
function FK_total = calculate_FK_force_by_component(active_list,nActive,...
    all_Area,all_Center,all_Tnorm,comp_start_idx,comp_end_idx,num_components,...
    elv,rho,g,amp,phase,kx,ky,dof,cg)

    FK_total = zeros(2,1);
    if isempty(kx) || isempty(amp), return; end

    for ia = 1:nActive
        comp = active_list(ia);
        if comp < 1 || comp > num_components, continue; end
        s = comp_start_idx(comp);
        e = comp_end_idx(comp);
        if s > e, continue; end

        area_comp = all_Area(s:e);
        center_comp = all_Center(s:e,:);
        center_comp(:,1) = center_comp(:,1) + dof(1) + cg(1);
        center_comp(:,2) = center_comp(:,2) + dof(2) + cg(2);
        tnorm_comp = all_Tnorm(s:e,:);
        elv_comp = elv(s:e);

        wp_comp = calcPressure(center_comp,elv_comp,amp,kx,ky,phase,rho,g);
        av_comp = tnorm_comp.*[area_comp area_comp];
        pressureVect = [wp_comp wp_comp].*(-av_comp);
        FK_total = FK_total + sum(pressureVect,1)';
    end
end

function wp = calcPressure(center,elv,amp,kx,ky,phase,rho,g)
    wp = zeros(size(center,1),1);
    z = center(:,3)-elv;
    for ic = 1:size(center,1)
        if z(ic) <= 0
            val = 0;
            for ik = 1:length(kx)
                k = sqrt(kx(ik)^2+ky(ik)^2);
                val = val + rho*g*amp(ik)* ...
                    cos(kx(ik)*center(ic,1)+ky(ik)*center(ic,2)+phase(ik))*exp(k*z(ic));
            end
            wp(ic) = val;
        end
    end
end

%% ========================================================================
% 子函数4：粘性阻力（逐部件）
% 使用每个部件瞬时 WetArea；波浪速度对所有波分量求和。
% 相对速度采用 U_rel = U_water - U_platform。
% ========================================================================
function Viscous_total = calculate_viscous_force_by_component(active_list,nActive,...
    TurbNum,all_Center,comp_start_idx,comp_end_idx,num_components,WetArea,vel,...
    rho,g,amp,phase,kx,ky,omega,Turb,dof,cg)

    Viscous_total = zeros(2,1);
    if isempty(kx) || isempty(amp), return; end

    for ia = 1:nActive
        comp = active_list(ia);
        if comp < 1 || comp > num_components, continue; end

        drag_coeff = 0.6;
        if comp <= length(Turb(TurbNum).ComponentInfo)
            drag_coeff = Turb(TurbNum).ComponentInfo(comp).DragCoeff;
        end

        proj_area = WetArea(comp);
        if proj_area <= 0, continue; end

        s = comp_start_idx(comp);
        e = comp_end_idx(comp);
        if s > e, continue; end
        center_comp = all_Center(s:e,:);
        center_comp(:,1) = center_comp(:,1) + dof(1) + cg(1);
        center_comp(:,2) = center_comp(:,2) + dof(2) + cg(2);
        % 用部件网格中心的平均 z 近似该部件的波浪粒子速度。
        avg_z = mean(center_comp(:,3));
        wave_vel = zeros(2,1);
        for ik = 1:length(kx)
            k = sqrt(kx(ik)^2+ky(ik)^2);
            if k > 0
                theta = kx(ik)*mean(center_comp(:,1)) + ...
                    ky(ik)*mean(center_comp(:,2)) + phase(ik);
                wave_amp_vel = amp(ik)*omega(ik)*exp(k*avg_z);
                wave_vel = wave_vel + wave_amp_vel*cos(theta)* ...
                    [kx(ik)/k;ky(ik)/k];
            end
        end

        % 流体相对于平台的 XY 相对速度
        % vel      = 平台速度
        % wave_vel = 波浪粒子速度
        % U_rel = U_water - U_platform
        rel_vel = wave_vel - vel;
        speed = norm(rel_vel);
        f_viscous = 0.5*rho*drag_coeff*proj_area*speed*rel_vel;

        % 保留原程序的数值保护，避免极端情况下把 Simulink 推爆。
        max_force = 1e6;
        f_norm = norm(f_viscous);
        if f_norm > max_force
            f_viscous = f_viscous/f_norm*max_force;
        end

        Viscous_total = Viscous_total + f_viscous;
    end
end

%% ========================================================================
% 子函数5：Radiation - 整体计算后按全部 WetArea 分配
% ========================================================================
function Rad_total = calculate_radiation_force_by_wet_area(...
    WetArea,active_list,nActive,vel_history,K_matrix,SampleTime)

    Rad_total = zeros(2,1);
    memory_term = zeros(2,1);

    if ~isempty(vel_history) && ~isempty(K_matrix)
        n_steps = min(size(vel_history,2),size(K_matrix,2));
        for i = 1:n_steps
            memory_term = memory_term + K_matrix(:,i).*vel_history(:,i)*SampleTime;
        end
    end
    Rad_force_total = -memory_term;

    total_wet_area = sum(WetArea);
    if total_wet_area <= 0, return; end

    for ia = 1:nActive
        comp = active_list(ia);
        if comp >= 1 && comp <= length(WetArea)
            ratio = WetArea(comp)/total_wet_area;
            Rad_total = Rad_total + Rad_force_total*ratio;
        end
    end
end

%% ========================================================================
% 子函数6：Diffraction - AQWA整体结果后按全部 WetArea 分配
% ========================================================================
function Diff_total = calculate_diffraction_force_by_wet_area(...
    wave_dat,diffData,dof,g,omega,WetArea,active_list,nActive)

    Diff_total = zeros(2,1);
    if isempty(wave_dat) || size(wave_dat,2) < 1, return; end
    if isempty(diffData.freqs), return; end

    Diff_total_6dof = calcDiffraction(wave_dat,diffData,dof,g,omega);
    if length(Diff_total_6dof) < 2, return; end

    total_wet_area = sum(WetArea);
    if total_wet_area <= 0, return; end

    for ia = 1:nActive
        comp = active_list(ia);
        if comp >= 1 && comp <= length(WetArea)
            ratio = WetArea(comp)/total_wet_area;
            Diff_total = Diff_total + Diff_total_6dof(1:2)*ratio;
        end
    end
end

%% ========================================================================
% 子函数7：AQWA Diffraction 整体计算
% ========================================================================
function F_diff_total = calcDiffraction(wave_dat,diffData,dof,g,omega)
    F_diff_total = zeros(6,1);
    [~,N_waves] = size(wave_dat);
    if N_waves < 1, return; end
    if isempty(diffData.freqs), return; end

    for j = 1:N_waves
        kx_j = wave_dat(1,j);
        ky_j = wave_dat(2,j);
        amp_j = wave_dat(3,j);
        phase_j = wave_dat(4,j);

        wave_angle_rad = atan2(ky_j,kx_j);
        wave_angle_deg = rad2deg(wave_angle_rad);

        [~,idx_f] = min(abs(diffData.freqs-omega(j)));
        [~,idx_d] = min(abs(diffData.dirs-wave_angle_deg));

        F_amp = squeeze(diffData.amp(idx_f,idx_d,:));
        F_phase_deg = squeeze(diffData.phase(idx_f,idx_d,:));
        F_phase_rad = deg2rad(F_phase_deg);

        phase_total = kx_j*dof(1) + ky_j*dof(2) + F_phase_rad + phase_j;
        F_j_raw = F_amp.*amp_j.*cos(phase_total);
        F_j_raw = F_j_raw(:);

        F_j = zeros(6,1);
        m = min(length(F_j_raw),6);
        for q = 1:m
            F_j(q) = F_j_raw(q);
        end
        F_diff_total = F_diff_total + F_j;
    end
end