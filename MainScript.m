repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end
cd(repoRoot);
addpath(genpath(repoRoot));
%% Setup workspace
clear;
clc;
close all;
format short g;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname','Times New Roman');
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);


%% Simulation settings
% OC4 Phase II load case 2.1:
% Regular waves | Support structure | No air | Regular airy: H = 6 m, T = 10 s
Case.Name = 'oc4_case_2p1_regular_wave';
Case.EnableWind = false;
Case.EnableWave = true;
Case.WaveType = 'regular_airy';
Case.NumTurbX = 1;
Case.NumTurbY = 1;
Case.WaveHeight = 11.12;  % m
Case.WavePeriod = 20.0;  % s
Case.WaveAmplitude = Case.WaveHeight/2;
Case.WaveDirectionDeg = 0;
Case.WavePhase0 = 0;
Case.SampleTime = 0.1;
% Case.TimeEnd = 3;  % Keep short for smoke testing; extend manually for paper comparisons.
Case.TimeEnd = 1000;
Case.WaveCaseDir = fullfile('.', 'e_oc4_case2p1');

% Simulation case
Global.TurbDia = 126;

if Case.EnableWind
    Global.MeanWindSpeed = 12;  % m/s
else
    Global.MeanWindSpeed = 0;
end
Farm.SpacingX = 7*Global.TurbDia;    %downstream spacing, m
Farm.SpacingY = 4*Global.TurbDia;


% Mesh settings
Solver.DnsElmLengthA = 1*Global.TurbDia;
Solver.NumMeshRedElm = 4;
Solver.CentralDiff = 1;

% Model dynamics?
Solver.TurbMotion = 1;

% Tuned parameters
Global.WakeExpConst = 0.08; 
Global.RefWindSpeed = max(Global.MeanWindSpeed, 1);
if Case.EnableWind
    Global.WakeExpRate = Global.WakeExpConst*Global.RefWindSpeed;
else
    Global.WakeExpRate = 0;
end
Global.GaussStd_Slope = 0.0250176061678923;
Global.GaussStd_Inter = 0.396350723166689;

% Actuator disk model
Solver.ActDiscModel = 'Vortex';
Solver.RotorMomModel = 'Jimenez';

% Fundamental sample time and impulse time setting
SampleTime = Case.SampleTime;
timeEnd = Case.TimeEnd;

%% Active component selector
% 唯一部件选择器：STL/AQWA 保持完整平台，只决定最终汇总哪些部件。
ActiveComponents = 1:19;


%%model selction
modelselect = 1;
%if modelselect equal to 1, the code would run as a 3 input model (yaw angle, blade pitch angle and generator torque)
% if the modelsect is not equal to 1, the code would run as a 3 input model (wind speed, yaw angle and AI factor)


%% Wind farm layout
% Wind farm layout with wind direction of 0 deg
Farm.NumTurbX = Case.NumTurbX;
Farm.NumTurbY = Case.NumTurbY;

% Total number of wind turbines
Farm.NumTurb = Farm.NumTurbX*Farm.NumTurbY;

%% Simulation inputs to vary
% Free stream wind conditions
InitFreeStreamWindVel = [Global.MeanWindSpeed 0]';
InitFreeStreamWindAccel = [0 0]';

%freeStreamImpulse = 0.45;
freeStreamImpulse = 0;
StepFreeStreamWindVel = [(Global.MeanWindSpeed + freeStreamImpulse) 0]';
StepFreeStreamWindAccel = [0 0]';

% Turbine inputs
InitAIFactorVec = 1/3*ones(Farm.NumTurb,1);
InitYawAngleVec = zeros(Farm.NumTurb,1);
InitRotorSpeedVec = zeros(Farm.NumTurb,1);
InitGeneTorqueVec = zeros(Farm.NumTurb,1);
InitBladePitch = zeros(Farm.NumTurb,1);

%warning
for i = 1:length(InitYawAngleVec)
    if (InitYawAngleVec(i)> 25*pi/180 || InitYawAngleVec(i)<-25*pi/180)
        disp('warning: values in initial yaw angle vector may need be converted to radians');
    end
end



%% External properties
% Gravitational acceleration
Global.GravAccel = 9.81;

% Fluid properties
Global.AirDensity = 1.225;
Global.WaterDensity = 1025;  % 文献值 1025 kg/m³

%% ========== 预定义重力加速度（用于波浪频率计算） ==========
Global.g = Global.GravAccel;

%% Floating wind turbine properties
% retardation function K的参数
k_length = 100; % irkbw.dat有存到多少个时间步之前的K矩阵
memory_length = 10; % 决定要让最近多少个时间步的记忆产生影响

for TurbNum = Farm.NumTurb:-1:1
    % Mechanical properties
    Turb(TurbNum).Mass = 13.444e6;  % 文献值 13,444,000 kg
    Turb(TurbNum).MomentInertia = 4.3703e7;
    
    % Turbine dimensions
    Turb(TurbNum).RotorDia = Global.TurbDia;
    Turb(TurbNum).RotorArea = pi/4*Turb(TurbNum).RotorDia^2;
    
    % Platform dimensions
    Turb(TurbNum).FLDistRelG  = 40.87;
    Turb(TurbNum).MidCylDia = 6.5;  %center cylinder
    Turb(TurbNum).TopCylDia = 12;
    Turb(TurbNum).BottCylDia = 24;
    Turb(TurbNum).MidCylLength = 20;
    Turb(TurbNum).TopCylLength = 14;
    Turb(TurbNum).BottCylLength = 6;
    
    % Platform hydrodynamic properties
    Turb(TurbNum).MidCylDragCoeff = 0.56;
    Turb(TurbNum).TopCylDragCoeff = 0.61;
    Turb(TurbNum).BottCylDragCoeff = 0.68;
    Turb(TurbNum).AddedMassCoeff = 0.63;
    Turb(TurbNum).EffDragFactor = 0.5*Global.WaterDensity*(...
        Turb(TurbNum).MidCylDragCoeff*Turb(TurbNum).MidCylDia*Turb(TurbNum).MidCylLength + ...
        3*Turb(TurbNum).TopCylDragCoeff*Turb(TurbNum).TopCylDia*Turb(TurbNum).TopCylLength + ...
        3*Turb(TurbNum).BottCylDragCoeff*Turb(TurbNum).BottCylDia*Turb(TurbNum).BottCylLength);
    Turb(TurbNum).AddedMass = pi/4*Turb(TurbNum).AddedMassCoeff*Global.WaterDensity*...
        (Turb(TurbNum).MidCylLength*Turb(TurbNum).MidCylDia^2 + ...
         3*Turb(TurbNum).TopCylLength*Turb(TurbNum).TopCylDia^2 + ...
         3*Turb(TurbNum).BottCylLength*Turb(TurbNum).BottCylDia^2);
    
    % 注意：水动力数据（area, center, tnorm, cg, K, diffData, component_info）
    % 将在后面的 "读取水动力数据" 部分统一读取
end

%% Wave propertiees
% =========================================================================
% 初始化与波浪数据预处理
% =========================================================================

% 波浪模式设置
% 0: 线性波 (Linear), 振幅固定，相位随时间线性变化
% 1: 非线性波 (Nonlinear), 每一步从文件读取振幅和相位
WaveDataMode = 0;

prepare_regular_wave_case(Case);
wave_case = [Case.WaveCaseDir filesep];

% 预读取数据 (无论哪种模式，这些基础数据只读一次)
disp('Loading Initial Wave Data...');

% --- 读取波数文件 (wave_number.dat) ---
fid_k = fopen([wave_case 'wave_number.dat'], 'r');
if fid_k == -1
    error('无法打开 wave_number.dat，请检查路径。');
end
data_k = fscanf(fid_k, '%f');
fclose(fid_k);

Wave.kx = data_k(1:2:end)';
Wave.ky = data_k(2:2:end)';

% --- 读取初始波面文件 (wave_0000.dat) ---
fid_0 = fopen([wave_case 'wave_0000.dat'], 'r');
if fid_0 == -1
    error('无法打开 wave_0000.dat，请检查路径。');
end
data_0 = fscanf(fid_0, '%f');
fclose(fid_0);

Wave.amp0 = data_0(1:2:end)';
Wave.phase0 = data_0(2:2:end)';

% --- 计算圆频率 omega ---
k_mag = sqrt(Wave.kx.^2 + Wave.ky.^2);
Wave.omega = sqrt(Global.g * k_mag);

Wave.wave_case = wave_case; 

disp(['Wave Data Loaded. Mode: ', num2str(WaveDataMode)]);
disp(['Number of Wave Components: ', num2str(length(Wave.kx))]);

%% ========== 读取水动力数据 ==========
% 您只需要修改下面两个路径
aqwa_path = 'E:\floating offshore wind turbine platforms\FOWFSimDyn1-wave_remix - test2.1-modular regular\read_AQWA\freq100(2)';
stl_path = 'E:\floating offshore wind turbine platforms\FOWFSimDyn1-wave_remix - test2.1-modular regular\read_AQWA\mesh';

% ============================================================
% 1. 读取全部 STL 部件
fprintf('\n========================================\n');
fprintf('步骤1: 读取全部 STL 部件（不按 ActiveComponents 删除）\n');
fprintf('========================================\n');
[component_mesh_data, part_names, all_vertices, all_faces] = load_all_stl_components(stl_path);
num_parts = length(component_mesh_data);
fprintf('STL 总部件数 = %d\n',num_parts);

ActiveComponents = unique(double(ActiveComponents(:)'),'stable');
if isempty(ActiveComponents), error('ActiveComponents 不能为空。'); end
for i=1:length(ActiveComponents)
    if ActiveComponents(i)<1 || ActiveComponents(i)>num_parts || ActiveComponents(i)~=floor(ActiveComponents(i))
        error('ActiveComponents 中存在非法部件编号：%g；STL 部件总数=%d。',ActiveComponents(i),num_parts);
    end
end
fprintf('ActiveComponents = ['); fprintf('%d ',ActiveComponents); fprintf(']\n');
fprintf('未选中部件仍保留在整体 AQWA/WetArea 基准中。\n');

for comp=1:num_parts
    if ~isscalar(component_mesh_data(comp).Volume)
        component_mesh_data(comp).Volume=sum(component_mesh_data(comp).Volume);
    end
end

% ============================================================
% 2. 读取 AQWA 水动力数据 (只读取频域系数)
% ============================================================
fprintf('\n========================================\n');
fprintf('步骤2: 读取AQWA水动力数据\n');
fprintf('========================================\n');

% 手动输入排水体积
displacement = 13986.8;  % m³, DeepCwind 文献 Table IX

for TurbNum = Farm.NumTurb:-1:1
    if TurbNum == Farm.NumTurb
        [dof, ~, ~, ~, ...
         Turb(TurbNum).cg, Global.rho, Global.g, Arrays.K, Turb(TurbNum).AddedMass, ...
         Turb(TurbNum).diffData, component_info] = ...
            read_hydro_data(SampleTime, k_length, stl_path, aqwa_path, displacement);

        fprintf('  AQWA Surge 方向附加质量 = %.2e kg\n',Turb(TurbNum).AddedMass);
        fprintf('  平台质量 = %.2e kg\n',Turb(TurbNum).Mass);
        fprintf('  总质量 = %.2e kg\n',Turb(TurbNum).Mass + Turb(TurbNum).AddedMass);
        fprintf('  排水体积 = %.1f m³ (手动输入)\n',displacement);

        fprintf('\n========================================\n');
        fprintf('步骤3: 构建部件信息（完整 STL，不使用 Enabled）\n');
        fprintf('========================================\n');
        component_info_simple=struct('ID',{},'Volume',{},'WetArea',{},'DragCoeff',{});
        for comp=1:num_parts
            component_info_simple(comp).ID=comp;
            component_info_simple(comp).Volume=double(component_mesh_data(comp).Volume);
            component_info_simple(comp).WetArea=sum(double(component_mesh_data(comp).Area));
            component_info_simple(comp).DragCoeff=0.6;
            fprintf('  部件 %d (%s): %d triangles, full area %.3f m^2, volume %.3f m^3\n', ...
                comp,component_mesh_data(comp).Name,length(component_mesh_data(comp).Area), ...
                component_info_simple(comp).WetArea,component_info_simple(comp).Volume);
        end

        Turb(TurbNum).ComponentInfo=component_info_simple;
        Turb(TurbNum).K=Arrays.K;
        Turb(TurbNum).VelHistory=zeros(2,k_length);
        Turb(TurbNum).VelHistoryPtr=1;
        Turb(TurbNum).ActiveComponents=ActiveComponents;
    else
        Turb(TurbNum).area=Turb(TurbNum+1).area;
        Turb(TurbNum).center=Turb(TurbNum+1).center;
        Turb(TurbNum).tnorm=Turb(TurbNum+1).tnorm;
        Turb(TurbNum).cg=Turb(TurbNum+1).cg;
        Turb(TurbNum).AddedMass=Turb(TurbNum+1).AddedMass;
        Turb(TurbNum).diffData=Turb(TurbNum+1).diffData;
        Turb(TurbNum).ComponentInfo=Turb(TurbNum+1).ComponentInfo;
        Turb(TurbNum).K=Turb(TurbNum+1).K;
        Turb(TurbNum).VelHistory=zeros(2,k_length);
        Turb(TurbNum).VelHistoryPtr=1;
        Turb(TurbNum).ActiveComponents=ActiveComponents;
    end
end

% ============================================================
% 4. 保存完整 STL 网格数据（含三角形三个顶点）
fprintf('\n保存 STLMeshData.mat（完整 STL + 顶点）...\n');
if exist('STLMeshData.mat','file'), delete('STLMeshData.mat'); end
num_parts=length(component_mesh_data);
total_mesh=0;
for comp=1:num_parts, total_mesh=total_mesh+length(component_mesh_data(comp).Area); end
all_Area=zeros(total_mesh,1,'double'); all_Center=zeros(total_mesh,3,'double'); all_Tnorm=zeros(total_mesh,2,'double'); all_Volume=zeros(num_parts,1,'double');
all_V1=zeros(total_mesh,3,'double'); all_V2=zeros(total_mesh,3,'double'); all_V3=zeros(total_mesh,3,'double');
comp_start_idx=zeros(num_parts,1,'double'); comp_end_idx=zeros(num_parts,1,'double');
idx_offset=0;
for comp=1:num_parts
    n_mesh=length(component_mesh_data(comp).Area); sidx=idx_offset+1; eidx=idx_offset+n_mesh;
    comp_start_idx(comp)=sidx; comp_end_idx(comp)=eidx;
    all_Area(sidx:eidx)=double(component_mesh_data(comp).Area);
    all_Center(sidx:eidx,:)=double(component_mesh_data(comp).Center);
    T=double(component_mesh_data(comp).Tnorm); all_Tnorm(sidx:eidx,:)=T(:,1:2);
    all_Volume(comp)=double(component_mesh_data(comp).Volume);
    Fblock=all_faces(sidx:eidx,:);
    all_V1(sidx:eidx,:)=all_vertices(Fblock(:,1),:);
    all_V2(sidx:eidx,:)=all_vertices(Fblock(:,2),:);
    all_V3(sidx:eidx,:)=all_vertices(Fblock(:,3),:);
    idx_offset=eidx;
end
save('STLMeshData.mat','all_Area','all_Center','all_Tnorm','all_Volume', ...
    'all_V1','all_V2','all_V3','comp_start_idx','comp_end_idx','num_parts','-v7.3');
fprintf('✓ 已保存 STLMeshData.mat: %d triangles, %d components\n',total_mesh,num_parts);

% 5. 初始化速度历史
% ============================================================
memory_length = size(Arrays.K, 2);
prev_vel_buffer = zeros(2, memory_length);
prev_buffer_ptr = 1;

fprintf('\n部件信息加载完成\n');
fprintf('  风机数量: %d\n', Farm.NumTurb);
fprintf('  每风机部件数: %d\n', length(Turb(1).ComponentInfo));

%% Mooring system properties
for TurbNum = Farm.NumTurb:-1:1
    % Cable mechanical properties
    Turb(TurbNum).CableNetDensity = 108.63;
    Turb(TurbNum).CableSpWeight = Global.GravAccel*Turb(TurbNum).CableNetDensity;
    Turb(TurbNum).CableStiff = 753.6e6;  
    Turb(TurbNum).SeabedFricCoeff = 1;
    
    % Cable dimensions
    Turb(TurbNum).CableLength = 835.5;
    Turb(TurbNum).CableDia = 0.0766;
    Turb(TurbNum).CableArea = pi*(Turb(TurbNum).CableDia^2)/4;
    
    % Mooring system dimensions
    Turb(TurbNum).VerDist = 186;
    Turb(TurbNum).AnchDistRelG = 837.6;
    Turb(TurbNum).HorDistInit = Turb(TurbNum).AnchDistRelG - Turb(TurbNum).FLDistRelG;
    Turb(TurbNum).AnchAngleRelX = 0;
end

%% Wind farm layout
TurbNum = Farm.NumTurb;
for TurbNumY = Farm.NumTurbY:-1:1
    for TurbNumX = Farm.NumTurbX:-1:1
        Turb(TurbNum).NeutralPosVec = zeros(2,1);
        Turb(TurbNum).NeutralPosVec(1) = 0 + (TurbNumX - 1)*Farm.SpacingX;
        Turb(TurbNum).NeutralPosVec(2) = 0 + (TurbNumY - 1)*Farm.SpacingY;
        TurbNum = TurbNum - 1;
    end
end 

Farm.AngleRelX = 0; % deg
Farm.CoordTransMat = [
    cosd(Farm.AngleRelX) -sind(Farm.AngleRelX);
    sind(Farm.AngleRelX) cosd(Farm.AngleRelX)];
TmpMat = zeros(2,Farm.NumTurb);
for TurbNum = 1:Farm.NumTurb
    Turb(TurbNum).NeutralPosVec = Farm.CoordTransMat*Turb(TurbNum).NeutralPosVec;
    TmpMat(:,TurbNum) = Turb(TurbNum).NeutralPosVec;
end
TmpMat = sortrows(TmpMat',1)';
for TurbNum = Farm.NumTurb:-1:1
    Turb(TurbNum).NeutralPosVec = TmpMat(:,TurbNum);
end

%% Calculate mooring system dimensions across wind farm
for TurbNum = Farm.NumTurb:-1:1
    Turb(TurbNum).AnchTransMat = [
        cosd(Turb(TurbNum).AnchAngleRelX) -sind(Turb(TurbNum).AnchAngleRelX);
        sind(Turb(TurbNum).AnchAngleRelX) cosd(Turb(TurbNum).AnchAngleRelX)];
    
    Turb(TurbNum).FL1PosRelG = Turb(TurbNum).FLDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(60);sind(60)];
    Turb(TurbNum).FL2PosRelG = Turb(TurbNum).FLDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(180);sind(180)];
    Turb(TurbNum).FL3PosRelG = Turb(TurbNum).FLDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(300);sind(300)];
    
    Turb(TurbNum).Anch1PosVec = Turb(TurbNum).AnchDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(60);sind(60)] + Turb(TurbNum).NeutralPosVec;
    Turb(TurbNum).Anch2PosVec = Turb(TurbNum).AnchDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(180);sind(180)] + Turb(TurbNum).NeutralPosVec;
    Turb(TurbNum).Anch3PosVec = Turb(TurbNum).AnchDistRelG*Farm.CoordTransMat*Turb(TurbNum).AnchTransMat*[cosd(300);sind(300)] + Turb(TurbNum).NeutralPosVec;
end

%% Generate lookup table for mooring line forces
Moor.HorDistSpacing = 1;

Moor.HorDistMin = 0;
Moor.HorDistCat = Turb(1).CableLength - Turb(1).VerDist;
Moor.HorForceTrans = (1 - (Turb(1).VerDist/Turb(1).CableLength - Turb(1).CableSpWeight*Turb(1).CableLength/(2*Turb(1).CableStiff))^2)*Turb(1).CableSpWeight*Turb(1).CableLength/(2*(Turb(1).VerDist/Turb(1).CableLength - Turb(1).CableSpWeight*Turb(1).CableLength/(2*Turb(1).CableStiff)));
Moor.VerForceTrans = Turb(1).CableSpWeight*Turb(1).CableLength;
Moor.HorDistTrans = (Moor.HorForceTrans/Turb(1).CableSpWeight)*(Turb(1).CableSpWeight*Turb(1).CableLength/Turb(1).CableStiff + asinh(Turb(1).CableSpWeight*Turb(1).CableLength/Moor.HorForceTrans));
Moor.MaxFac = 1.5;
Moor.HorDistMax = Moor.MaxFac*sqrt(Turb(1).CableLength^2 - Turb(1).VerDist^2);
HorDistVec = (0:Moor.HorDistSpacing:Moor.HorDistMax)';

HorForceVec = zeros(length(HorDistVec),1);
VerForceVec = zeros(length(HorDistVec),1);
Moor.HorForceGuess = Moor.HorForceTrans;
Moor.VerForceGuess = Moor.VerForceTrans;
Moor_fsolveOptions = optimoptions('fsolve',...
    'Display','off',...
    'FiniteDifferenceType','central',...
    'StepTolerance',1e-16,...
    'FunctionTolerance',1e-16,...
    'MaxFunctionEvaluations',10000,...
    'MaxIterations',10000,...
    'OptimalityTolerance',1e-16);
for i = 1:length(HorDistVec)
    if HorDistVec(i) <= Moor.HorDistCat
        HorForceVec(i) = 0;
        VerForceVec(i) = Turb(1).CableSpWeight*Turb(1).VerDist;
    elseif HorDistVec(i) <= Moor.HorDistTrans
        Func = @(x) [
            (-HorDistVec(i) + (Turb(1).CableLength - x(2)/Turb(1).CableSpWeight +((1 + x(1)/Turb(1).CableStiff)^3 - ((1 + x(1)/Turb(1).CableStiff)^2 - 2*Turb(1).CableSpWeight*Turb(1).SeabedFricCoeff/Turb(1).CableStiff*min(Turb(1).CableLength - x(2)/Turb(1).CableSpWeight,x(1)*(1 + x(1)/(2*Turb(1).CableStiff))/(Turb(1).CableSpWeight*Turb(1).SeabedFricCoeff)))^(3/2))/(3*Turb(1).SeabedFricCoeff*Turb(1).CableSpWeight/Turb(1).CableStiff) - min(Turb(1).CableLength - x(2)/Turb(1).CableSpWeight,x(1)*(1 + x(1)/(2*Turb(1).CableStiff))/(Turb(1).CableSpWeight*Turb(1).SeabedFricCoeff))) + x(1)/Turb(1).CableSpWeight*(x(2)/Turb(1).CableStiff + asinh(x(2)/x(1))));
            (-Turb(1).VerDist + (1/Turb(1).CableSpWeight)*((x(2)^2)/(2*Turb(1).CableStiff) - x(1)*(1 - sqrt(1 + (x(2)/x(1))^2))))];
        Vars = fsolve(Func,[Moor.HorForceGuess Moor.VerForceGuess]',Moor_fsolveOptions);
        HorForceVec(i) = Vars(1);
        VerForceVec(i) = Vars(2);
    elseif HorDistVec(i) > Moor.HorDistTrans
        Func = @(x) [
            (-HorDistVec(i) + x(1)/Turb(1).CableSpWeight*(Turb(1).CableSpWeight*Turb(1).CableLength/Turb(1).CableStiff + asinh(x(2)/x(1)) - asinh((x(2) - Turb(1).CableSpWeight*Turb(1).CableLength)/x(1))));
            (-Turb(1).VerDist + (Turb(1).CableLength/Turb(1).CableStiff*(x(2) - Turb(1).CableSpWeight*...
            Turb(1).CableLength/2) + x(1)/Turb(1).CableSpWeight*(sqrt(1 + (x(2)/x(1))^2) - sqrt(1 + ((x(2) - Turb(1).CableSpWeight*Turb(1).CableLength)/x(1))^2))))];
        Vars = fsolve(Func,[Moor.HorForceGuess Moor.VerForceGuess]',Moor_fsolveOptions);
        HorForceVec(i) = Vars(1);
        VerForceVec(i) = Vars(2);
    end
end

HorDistData = HorDistVec;
HorForceData = HorForceVec;
VerForceData = VerForceVec;

%% Build finite difference mesh
Farm.ExtraSpacingX = 7*Global.TurbDia;
Farm.ExtraSpacingY = 3*Global.TurbDia;
Farm.MaxNumWakePoints = 1000;
Farm.LengthX = max(TmpMat(1,:)) - min(TmpMat(1,:));
Farm.LengthY = max(TmpMat(2,:)) - min(TmpMat(2,:));
Farm.TotalNumWakePoints = 0;
for TurbNum = Farm.NumTurb:-1:1
    Turb(TurbNum).DnsLength = Farm.LengthX*2 + Farm.ExtraSpacingX - (Turb(TurbNum).NeutralPosVec(1) - Turb(1).NeutralPosVec(1));
    Turb(TurbNum).NumElm_A = ceil(Turb(TurbNum).DnsLength/Solver.DnsElmLengthA);
    Turb(TurbNum).NumElm = Turb(TurbNum).NumElm_A + Solver.NumMeshRedElm;
    Turb(TurbNum).ElmLength_A = Turb(TurbNum).DnsLength/Turb(TurbNum).NumElm_A;
    Turb(TurbNum).NumWakePoints = Turb(TurbNum).NumElm;
    Farm.TotalNumWakePoints = Farm.TotalNumWakePoints + Turb(TurbNum).NumWakePoints;
    Turb(TurbNum).WakePointPos_X = zeros(1,Farm.MaxNumWakePoints);
    for PointNum = 1:Turb(TurbNum).NumElm_A
        Turb(TurbNum).WakePointPos_X(PointNum) = PointNum*Turb(TurbNum).ElmLength_A;
    end
    ElmLengthTmp = Turb(TurbNum).ElmLength_A/2;
    for PointNum = Turb(TurbNum).NumElm_A + 1:Turb(TurbNum).NumWakePoints
        Turb(TurbNum).WakePointPos_X(PointNum) = Turb(TurbNum).WakePointPos_X(PointNum - 1) + ElmLengthTmp;
        ElmLengthTmp = ElmLengthTmp/2;
    end
end

%% Initialize dynamic arrays
Global.FreeStreamWindVel = InitFreeStreamWindVel;
Global.FreeStreamWindAccel = InitFreeStreamWindAccel;
Global.FreeStreamUnitVec = zeros(2,1);
for TurbNum = Farm.NumTurb:-1:1
    Turb(TurbNum).PosVec = zeros(2,1);
    Turb(TurbNum).VelVec = zeros(2,1);
    Turb(TurbNum).WakePointPos_Y = zeros(1,Farm.MaxNumWakePoints);
    Turb(TurbNum).WakePointVelVec = zeros(2,Farm.MaxNumWakePoints);
    Turb(TurbNum).WakePointDia = zeros(1,Farm.MaxNumWakePoints);
    Turb(TurbNum).AIFactor = 0;
    Turb(TurbNum).YawAngle = 0;
    Turb(TurbNum).Power = 0;
    Turb(TurbNum).ThrustCoeff = 0;
    Turb(TurbNum).PowerCoeff = 0;
    Turb(TurbNum).RotorNormVec = zeros(2,1);
    Turb(TurbNum).ThrustForceVec = zeros(2,1);
    Turb(TurbNum).UpsWindVelVec = zeros(2,1);
    Turb(TurbNum).RelUpsWindVel = zeros(2,1);
    Turb(TurbNum).HydroForceVec = zeros(2,1);
    Turb(TurbNum).MoorForceVec = zeros(2,1);
    Turb(TurbNum).ForceVec = zeros(2,1);
    Turb(TurbNum).AccelVec = zeros(2,1);
    Turb(TurbNum).RelWindAngle = 0;
    Turb(TurbNum).RelYawAngle = 0;
    Turb(TurbNum).InitWakeVelVec = zeros(2,1);
    Turb(TurbNum).RelUpsWindVelVec = zeros(2,1);
    Turb(TurbNum).BladePitch=0;
    Turb(TurbNum).GeneTorque=0;
    Turb(TurbNum).RotorSpeed=0;
    Turb(TurbNum).NormalWindVel=0;
    Turb(TurbNum).Lambda=0;
    Turb(TurbNum).TurbTorque=0;
end
Arrays.TurbPosStateDerivMat = zeros(2,Farm.NumTurb);
Arrays.TurbVelStateDerivMat = zeros(2,Farm.NumTurb);
Arrays.WakeStateDerivMat = zeros(4,Farm.TotalNumWakePoints);
Arrays.TurbPowerVec = zeros(1,Farm.NumTurb);
Arrays.WakeLocalVel = zeros(1,Farm.MaxNumWakePoints);
Arrays.WakeLocalAccel = zeros(2,Farm.MaxNumWakePoints);
Arrays.WakeLocalDiaROT = zeros(1,Farm.MaxNumWakePoints);
Arrays.Output_PosX = zeros(Farm.NumTurb,1);
Arrays.Output_PosY = zeros(Farm.NumTurb,1);
Arrays.Output_Power = zeros(Farm.NumTurb,1);



%% Initial conditions
InitTurbPosMat = zeros(2,Farm.NumTurb);
InitTurbVelMat = zeros(2,Farm.NumTurb);
InitWakeStateMat = zeros(4,Farm.TotalNumWakePoints);

BaseNum = 0;
for TurbNum = 1:Farm.NumTurb
    InitTurbPosMat(:,TurbNum) = Turb(TurbNum).NeutralPosVec;
    InitTurbVelMat(:,TurbNum) = [0 0]';
    
    if Case.EnableWind && norm(InitFreeStreamWindVel) > 0
        InitWakeStateMat(:,BaseNum + 1:BaseNum + Turb(TurbNum).NumWakePoints) = [
            [0 1]*InitFreeStreamWindVel/([1 0]*InitFreeStreamWindVel)*Turb(TurbNum).WakePointPos_X(1:Turb(TurbNum).NumWakePoints);
            (InitFreeStreamWindVel - InitTurbVelMat(:,TurbNum))*ones(1,Turb(TurbNum).NumWakePoints);
            Turb(TurbNum).RotorDia + Global.WakeExpRate/norm(InitFreeStreamWindVel)*Turb(TurbNum).WakePointPos_X(1:Turb(TurbNum).NumWakePoints)];
    else
        InitWakeStateMat(:,BaseNum + 1:BaseNum + Turb(TurbNum).NumWakePoints) = [
            zeros(1,Turb(TurbNum).NumWakePoints);
            zeros(2,Turb(TurbNum).NumWakePoints);
            Turb(TurbNum).RotorDia*ones(1,Turb(TurbNum).NumWakePoints)];
    end
    
    BaseNum = BaseNum + Turb(TurbNum).NumWakePoints;
end

InitStateVec = [
    reshape(InitTurbPosMat,[],1);
    reshape(InitTurbVelMat,[],1);
    reshape(InitWakeStateMat,[],1)];
StateDerivs = zeros(size(InitStateVec));

%% Run simulation
% ============================================================
% 保存 Parameters_Simu.mat（不含 Enabled）
fprintf('保存 Parameters_Simu.mat ...\n');
for TurbNum=1:Farm.NumTurb
    nComp=length(Turb(TurbNum).ComponentInfo);
    clean_info=struct('ID',{},'Volume',{},'WetArea',{},'DragCoeff',{});
    for comp=1:nComp
        clean_info(comp).ID=Turb(TurbNum).ComponentInfo(comp).ID;
        clean_info(comp).Volume=Turb(TurbNum).ComponentInfo(comp).Volume;
        clean_info(comp).WetArea=Turb(TurbNum).ComponentInfo(comp).WetArea;
        clean_info(comp).DragCoeff=Turb(TurbNum).ComponentInfo(comp).DragCoeff;
    end
    Turb(TurbNum).ComponentInfo=clean_info;
    Turb(TurbNum).ActiveComponents=ActiveComponents;
end
WaveData.kx=Wave.kx; WaveData.ky=Wave.ky; WaveData.amp0=Wave.amp0; WaveData.phase0=Wave.phase0; WaveData.omega=Wave.omega; WaveData.wave_case=Wave.wave_case;
save('Parameters_Simu.mat','Solver','Global','Farm','Turb','Moor','Arrays','Case','WaveData','ActiveComponents');
fprintf('✓ 已保存 Parameters_Simu.mat\n');
fprintf('ActiveComponents = ['); fprintf('%d ',ActiveComponents); fprintf(']\n');

InitFreeStreamWindVel = Global.MeanWindSpeed;
StepFreeStreamWindVel = Global.MeanWindSpeed + freeStreamImpulse;

% 运行仿真
clear test_cal_FK_diff;
simOut = sim('Simulink_StateDerivVec_v8_ViscousSignFixed');

%% Plot
time = simOut.tout;

PosX = simOut.PosX.signals.values;
PosY = simOut.PosY.signals.values;
Power = simOut.Power.signals.values;

figure(1)

subplot(3,1,1)
for i = 1:Farm.NumTurb
    stairs(time, PosX(:,i)); hold on;
end
title('Position X');
xlabel('Time [s]');
ylabel('Position [m]');
legend('Turbine1');
grid on


subplot(3,1,2)
for i = 1:Farm.NumTurb
    stairs(time, PosY(:,i)); hold on;
end
title('Position Y');
xlabel('Time [s]');
ylabel('Position [m]');
Legend = cell(Farm.NumTurb,1);
for iter=1:Farm.NumTurb
    Legend{iter}=strcat('Turbine', num2str(iter));
end
legend(Legend);
grid on

subplot(3,1,3)
for i = 1:Farm.NumTurb
    stairs(time, Power(:,i)); hold on;
end
title('Power');
xlabel('Time [s]');
ylabel('Power [W]');
Legend = cell(Farm.NumTurb,1);
for iter=1:Farm.NumTurb
    Legend{iter}=strcat('Turbine', num2str(iter));
end
legend(Legend);
grid on

%% ============================================================
% 计算并输出 RAO 值
% ============================================================
% 获取波幅
if isfield(Case, 'WaveAmplitude') && Case.WaveAmplitude > 0
    wave_amp = Case.WaveAmplitude;
else
    wave_amp = Case.WaveHeight / 2;
end

% 取稳定段（后 50% 数据）
total_len = length(time);
steady_start_idx = round(total_len * 0.5);
steady_PosX = PosX(steady_start_idx:end, 1);
steady_PosY = PosY(steady_start_idx:end, 1);

% 计算响应幅值（峰峰值的一半）
surge_amp = (max(steady_PosX) - min(steady_PosX)) / 2;
sway_amp = (max(steady_PosY) - min(steady_PosY)) / 2;

% 计算 RAO
RAO_surge = surge_amp / wave_amp;
RAO_sway = sway_amp / wave_amp;

% 输出 RAO 值
fprintf('\n========================================\n');
fprintf('RAO 计算结果:\n');
fprintf('========================================\n');
fprintf('  波幅 A = %.4f m\n', wave_amp);
if isfield(Case, 'WavePeriod')
    fprintf('  波周期 T = %.2f s\n', Case.WavePeriod);
end
fprintf('  Surge 响应幅值 = %.4f m\n', surge_amp);
fprintf('  Surge RAO = %.4f (m/m)\n', RAO_surge);
fprintf('  Sway 响应幅值 = %.4f m\n', sway_amp);
fprintf('  Sway RAO = %.4f (m/m)\n', RAO_sway);
fprintf('========================================\n');

%% ========================================================================
% Local function: 读取全部 STL 部件
% ========================================================================
function [component_mesh_data,part_names,all_vertices,all_faces] = load_all_stl_components(stl_path)
files=dir(fullfile(stl_path,'*.stl'));
if isempty(files), error('STL 文件夹中没有找到 .stl 文件：%s',stl_path); end
names={files.name}; [~,order]=sort(lower(names)); files=files(order);
num_parts=length(files); part_names=cell(num_parts,1); all_vertices=zeros(0,3); all_faces=zeros(0,3); offsetV=0; part_faces=cell(num_parts,1);
for comp=1:num_parts
    part_names{comp}=files(comp).name; fn=fullfile(stl_path,files(comp).name);
    try
        TR=stlread(fn); V=double(TR.Points); F=double(TR.ConnectivityList);
    catch ME
        error('读取 STL 失败：Part %d (%s): %s',comp,fn,ME.message);
    end
    part_faces{comp}=F; all_vertices=[all_vertices;V]; all_faces=[all_faces;F+offsetV]; offsetV=offsetV+size(V,1);
end
total_tri=size(all_faces,1); area=zeros(total_tri,1); normal=zeros(total_tri,3); center=zeros(total_tri,3); offset=0;
for comp=1:num_parts
    n=size(part_faces{comp},1);
    for j=1:n
        ii=offset+j; v1=all_vertices(all_faces(ii,1),:); v2=all_vertices(all_faces(ii,2),:); v3=all_vertices(all_faces(ii,3),:);
        nv=cross(v2-v1,v3-v1); nm=norm(nv);
        if nm>0, normal(ii,:)=nv/nm; area(ii)=0.5*nm; else, normal(ii,:)=[0 0 1]; area(ii)=0; end
        center(ii,:)=(v1+v2+v3)/3;
    end
    offset=offset+n;
end
component_mesh_data=struct('ID',{},'Name',{},'Area',{},'Center',{},'Tnorm',{},'Volume',{}); offset=0;
for comp=1:num_parts
    n=size(part_faces{comp},1); idx=offset+(1:n);
    component_mesh_data(comp).ID=comp; component_mesh_data(comp).Name=part_names{comp}; component_mesh_data(comp).Area=area(idx); component_mesh_data(comp).Center=center(idx,:); component_mesh_data(comp).Tnorm=normal(idx,:);
    Vsum=0;
    for j=1:n
        ii=offset+j; v1=all_vertices(all_faces(ii,1),:); v2=all_vertices(all_faces(ii,2),:); v3=all_vertices(all_faces(ii,3),:); Vsum=Vsum+dot(cross(v1,v2),v3)/6;
    end
    component_mesh_data(comp).Volume=abs(Vsum); offset=offset+n;
end
end
