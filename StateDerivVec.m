function [StateDerivs,Output_PosX,Output_PosY,Output_Power,Output_PowerC,Output_ThrustC,Output_RotorAcc,Output_RotorSpeedv,shutdown,windx,new_vel_buffer,new_buffer_ptr] = ...
    StateDerivVec(WaveExcitationForce, SampleTime, modelselect, HorDistData, HorForceData, Time, StateVec, ...
    AIFactorVec, YawAngleVec, FreeStreamWindVel, GeneTorqueVec, BladepitchVec, RotorSpeedVec, ...
    dRotorSpeedVec, preshut, prev_vel_buffer, prev_buffer_ptr)

    persistent last_hydro_diag_time
    if isempty(last_hydro_diag_time)
        last_hydro_diag_time = -inf;
    end

    %% Load structures
    Pars = load('Parameters_Simu.mat');

    Solver = Pars.Solver;
    Global = Pars.Global;
    Farm = Pars.Farm;
    Turb = Pars.Turb;
    Arrays = Pars.Arrays;
    TurbPosStateDerivMat = Arrays.TurbPosStateDerivMat;
    TurbVelStateDerivMat = Arrays.TurbVelStateDerivMat;
    WakeStateDerivMat = Arrays.WakeStateDerivMat;
    TurbPowerVec = Arrays.TurbPowerVec;
    WakeLocalVel = Arrays.WakeLocalVel;
    WakeLocalAccel = Arrays.WakeLocalAccel;
    WakeLocalDiaROT = Arrays.WakeLocalDiaROT;
    Output_PosX = Arrays.Output_PosX;
    Output_PosY = Arrays.Output_PosY;
    Output_Power = Arrays.Output_Power;
    Output_ThrustC = Arrays.Output_Power;
    Output_PowerC = Arrays.Output_Power;
    Output_RotorAcc = Arrays.Output_Power;
    Output_RotorSpeedv = Arrays.Output_Power;
    shutdown = Arrays.Output_Power;
    windx = Arrays.Output_Power;

    enableWind = false;
    if isfield(Pars, 'Case') && isfield(Pars.Case, 'EnableWind')
        enableWind = Pars.Case.EnableWind;
    end

    memory_length = size(prev_vel_buffer, 2);
    if memory_length == 0
        memory_length = 10;
        prev_vel_buffer = zeros(2, memory_length);
    end
    new_vel_buffer = prev_vel_buffer;
    new_buffer_ptr = prev_buffer_ptr;

    data = load('paraVal.mat');
    Lambda.mesh = data.paraVal.aero.lambda_grid;
    Beta.mesh = data.paraVal.aero.beta_grid;
    Cp.mesh = data.paraVal.aero.cp_grid;
    Ct.mesh = data.paraVal.aero.ct_grid;

    %% Extract states and inputs
    TurbPosMat = reshape(StateVec(1:2*Farm.NumTurb), 2, Farm.NumTurb);
    TurbVelMat = reshape(StateVec(2*Farm.NumTurb + 1:4*Farm.NumTurb), 2, Farm.NumTurb);
    WakeStateMat = reshape(StateVec(4*Farm.NumTurb + 1:4*Farm.NumTurb + 4*Farm.TotalNumWakePoints), 4, Farm.TotalNumWakePoints);
    for TurbNum = 1:Farm.NumTurb
        Turb(TurbNum).PosVec = TurbPosMat(:, TurbNum);
        Turb(TurbNum).VelVec = TurbVelMat(:, TurbNum);

        Turb(TurbNum).WakePointPos_Y(1:Turb(TurbNum).NumWakePoints) = WakeStateMat(1, 1:Turb(TurbNum).NumWakePoints);
        Turb(TurbNum).WakePointVelVec(:, 1:Turb(TurbNum).NumWakePoints) = WakeStateMat(2:3, 1:Turb(TurbNum).NumWakePoints);
        Turb(TurbNum).WakePointDia(1:Turb(TurbNum).NumWakePoints) = WakeStateMat(4, 1:Turb(TurbNum).NumWakePoints);
        WakeStateMat(:, 1:Turb(TurbNum).NumWakePoints) = [];

        Turb(TurbNum).YawAngle = YawAngleVec(TurbNum);
        if modelselect == 1
            Turb(TurbNum).BladePitch = BladepitchVec(TurbNum);
            Turb(TurbNum).GeneTorque = GeneTorqueVec(TurbNum);
            if Time == 0
                Turb(TurbNum).RotorSpeed = RotorSpeedVec(TurbNum);
            else
                Turb(TurbNum).RotorSpeed = RotorSpeedVec(TurbNum) + dRotorSpeedVec(TurbNum);
            end
        else
            Turb(TurbNum).AIFactor = AIFactorVec(TurbNum);
        end
    end

    % Free stream wind condition
    Global.FreeStreamWindVel = FreeStreamWindVel;
    Global.FreeStreamWindAccel = [0 0]';

    %% Define free stream wind unit vector
    if norm(Global.FreeStreamWindVel) > 0
        Global.FreeStreamUnitVec = Global.FreeStreamWindVel / norm(Global.FreeStreamWindVel);
    else
        Global.FreeStreamUnitVec = [1; 0];
    end

    %% Begin loop for state derivatives
    BaseNum = 0;
    for TurbNum = 1:Farm.NumTurb

        %% Calculate turbine platform state-derivatives
        if enableWind
            if TurbNum == 1
                Turb(TurbNum).UpsWindVelVec = Global.FreeStreamWindVel;
            else
                EffVelDeficit = 0;
                for UpsTurbNum = 1:TurbNum - 1
                    if preshut(TurbNum - 1) == 1 && UpsTurbNum > 1
                        UpsTurbNum1 = UpsTurbNum - 1;
                    else
                        UpsTurbNum1 = UpsTurbNum;
                    end

                    TurbPosXRelUpsTurb = Turb(TurbNum).PosVec(1) - Turb(UpsTurbNum1).PosVec(1);
                    if TurbPosXRelUpsTurb >= 0
                        UpsWakePosXData = [0 Turb(UpsTurbNum1).WakePointPos_X(1:Turb(TurbNum).NumWakePoints)];
                        UpsWakePosYData = [0 Turb(UpsTurbNum1).WakePointPos_Y(1:Turb(TurbNum).NumWakePoints)];
                        UpsWakeVelXData = [Turb(UpsTurbNum1).InitWakeVelVec(1) Turb(UpsTurbNum1).WakePointVelVec(1, 1:Turb(TurbNum).NumWakePoints)];
                        UpsWakeVelYData = [Turb(UpsTurbNum1).InitWakeVelVec(2) Turb(UpsTurbNum1).WakePointVelVec(2, 1:Turb(TurbNum).NumWakePoints)];
                        UpsWakeDiaData = [Turb(UpsTurbNum1).RotorDia Turb(UpsTurbNum1).WakePointDia(1:Turb(TurbNum).NumWakePoints)];

                        for index = 1:length(UpsWakeVelXData)
                            if isnan(UpsWakeVelXData(index)) || isinf(UpsWakeVelXData(index))
                                UpsWakeVelXData(index) = 0;
                            end
                        end
                        for indexs = 1:length(UpsWakeVelYData)
                            if isnan(UpsWakeVelYData(indexs)) || isinf(UpsWakeVelYData(indexs))
                                UpsWakeVelYData(indexs) = 0;
                            end
                        end

                        UpsWakePosY = Turb(UpsTurbNum1).PosVec(2) + interp1(UpsWakePosXData, UpsWakePosYData, TurbPosXRelUpsTurb, 'linear', 0);
                        UpsWakeVelX = Turb(UpsTurbNum1).VelVec(1) + interp1(UpsWakePosXData, UpsWakeVelXData, TurbPosXRelUpsTurb, 'linear', 0);
                        UpsWakeVelY = Turb(UpsTurbNum1).VelVec(2) + interp1(UpsWakePosXData, UpsWakeVelYData, TurbPosXRelUpsTurb, 'linear', 0);
                        UpsWakeDia = interp1(UpsWakePosXData, UpsWakeDiaData, TurbPosXRelUpsTurb, 'linear', 0);
                        UpsWakeOffset = abs(UpsWakePosY - Turb(TurbNum).PosVec(2));
                        if UpsWakeOffset >= Turb(TurbNum).RotorDia/2 + UpsWakeDia/2
                            EffVelVec = Global.FreeStreamWindVel;
                        else
                            b = Turb(TurbNum).RotorDia*(Global.GaussStd_Slope*TurbPosXRelUpsTurb/Turb(TurbNum).RotorDia + Global.GaussStd_Inter);
                            a_x = 0.5*(((UpsWakeDia/2)/b)^2)*([1 0]*Global.FreeStreamWindVel - UpsWakeVelX);
                            a_y = 0.5*(((UpsWakeDia/2)/b)^2)*([0 1]*Global.FreeStreamWindVel - UpsWakeVelY);
                            UpsWakeVelX_Gauss = @(r) [1 0]*Global.FreeStreamWindVel - a_x*exp(-(r.^2)/(2*(b^2)));
                            UpsWakeVelY_Gauss = @(r) [0 1]*Global.FreeStreamWindVel - a_y*exp(-(r.^2)/(2*(b^2)));
                            Func_x = @(r) r.*AngleOLFunc(UpsWakeOffset, r, Turb(TurbNum).RotorDia/2).*UpsWakeVelX_Gauss(r);
                            Func_y = @(r) r.*AngleOLFunc(UpsWakeOffset, r, Turb(TurbNum).RotorDia/2).*UpsWakeVelY_Gauss(r);
                            r_vec = linspace(max(UpsWakeOffset - Turb(UpsTurbNum1).RotorDia/2, 0), UpsWakeOffset + Turb(TurbNum).RotorDia/2, 10);
                            EffVelVec = [
                                trapz(r_vec, Func_x(r_vec))/Turb(TurbNum).RotorArea;
                                trapz(r_vec, Func_y(r_vec))/Turb(TurbNum).RotorArea];
                        end
                        EffVelDeficit = EffVelDeficit + (norm(Global.FreeStreamWindVel) - dot(EffVelVec, Global.FreeStreamUnitVec))^2;
                    end
                end
                Turb(TurbNum).UpsWindVelVec = (norm(Global.FreeStreamWindVel) - sqrt(EffVelDeficit))*Global.FreeStreamUnitVec;
                if Turb(TurbNum).UpsWindVelVec(1) < 0
                    Turb(TurbNum).UpsWindVelVec(1) = 0;
                end
                if Turb(TurbNum).UpsWindVelVec(2) < 0
                    Turb(TurbNum).UpsWindVelVec(2) = 0;
                end
            end
        else
            Turb(TurbNum).UpsWindVelVec = [0; 0];
        end

        %% Ali's model to find Power and Thrust coefficients
        if enableWind
            Turb(TurbNum).RelUpsWindVelVec = Turb(TurbNum).UpsWindVelVec - Turb(TurbNum).VelVec;
            Turb(TurbNum).RelWindAngle = atan(Turb(TurbNum).RelUpsWindVelVec(2)/Turb(TurbNum).RelUpsWindVelVec(1));
            Turb(TurbNum).RelYawAngle = Turb(TurbNum).YawAngle - Turb(TurbNum).RelWindAngle;

            if modelselect == 1
                Turb(TurbNum).NormalWindVel = norm(Turb(TurbNum).RelUpsWindVelVec)*cos(Turb(TurbNum).RelYawAngle);
                Turb(TurbNum).Lambda = Turb(TurbNum).RotorSpeed*Turb(TurbNum).RotorDia/2/Turb(TurbNum).NormalWindVel;
                Turb(TurbNum).PowerCoeff = interp2(Beta.mesh, Lambda.mesh, Cp.mesh, Turb(TurbNum).BladePitch, Turb(TurbNum).Lambda);
                Turb(TurbNum).ThrustCoeff = interp2(Beta.mesh, Lambda.mesh, Ct.mesh, Turb(TurbNum).BladePitch, Turb(TurbNum).Lambda);
                Turb(TurbNum).TurbTorque = 0.5*Global.AirDensity*Turb(TurbNum).RotorArea*(Turb(TurbNum).NormalWindVel)^3*Turb(TurbNum).PowerCoeff/Turb(TurbNum).RotorSpeed;
                Output_RotorAcc(TurbNum) = (Turb(TurbNum).TurbTorque - Turb(TurbNum).GeneTorque)/Turb(TurbNum).MomentInertia;
                if isinf(Output_RotorAcc(TurbNum)) || isnan(Output_RotorAcc(TurbNum))
                    Output_RotorAcc(TurbNum) = 0;
                end
            else
                if Turb(TurbNum).UpsWindVelVec(1) <= 3
                    Turb(TurbNum).AIFactor = 0;
                else
                    Turb(TurbNum).AIFactor = AIFactorVec(TurbNum);
                end
                if strcmp(Solver.ActDiscModel, 'Momentum')
                    Turb(TurbNum).ThrustCoeff = 4*Turb(TurbNum).AIFactor*(cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor);
                    Turb(TurbNum).PowerCoeff = 4*Turb(TurbNum).AIFactor*((cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor)^2);
                elseif strcmp(Solver.ActDiscModel, 'Glauert')
                    Turb(TurbNum).ThrustCoeff = 4*Turb(TurbNum).AIFactor*sqrt(1 - Turb(TurbNum).AIFactor*(2*cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor));
                    Turb(TurbNum).PowerCoeff = 4*Turb(TurbNum).AIFactor*(cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor)*sqrt(1 - Turb(TurbNum).AIFactor*(2*cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor));
                elseif strcmp(Solver.ActDiscModel, 'Vortex')
                    DefAngle = (0.6*Turb(TurbNum).AIFactor + 1)*Turb(TurbNum).RelYawAngle;
                    Turb(TurbNum).ThrustCoeff = 4*Turb(TurbNum).AIFactor*(cos(Turb(TurbNum).RelYawAngle) + tan(DefAngle/2)*sin(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor*((sec(DefAngle/2))^2));
                    Turb(TurbNum).PowerCoeff = 4*Turb(TurbNum).AIFactor*(cos(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor)*(cos(Turb(TurbNum).RelYawAngle) + tan(DefAngle/2)*sin(Turb(TurbNum).RelYawAngle) - Turb(TurbNum).AIFactor*((sec(DefAngle/2))^2));
                else
                    error('Selected actuator disc model does not exist...');
                end
            end

            Turb(TurbNum).RotorNormVec = [
                cos(Turb(TurbNum).YawAngle);
                sin(Turb(TurbNum).YawAngle)];
            Turb(TurbNum).ThrustForceVec = ...
                0.5*Turb(TurbNum).ThrustCoeff*Global.AirDensity*Turb(TurbNum).RotorArea* ...
                (norm(Turb(TurbNum).RelUpsWindVelVec)^2)*Turb(TurbNum).RotorNormVec;
        else
            Turb(TurbNum).RelUpsWindVelVec = [0; 0];
            Turb(TurbNum).RelWindAngle = 0;
            Turb(TurbNum).RelYawAngle = Turb(TurbNum).YawAngle;
            Turb(TurbNum).NormalWindVel = 0;
            Turb(TurbNum).Lambda = 0;
            Turb(TurbNum).PowerCoeff = 0;
            Turb(TurbNum).ThrustCoeff = 0;
            Turb(TurbNum).TurbTorque = 0;
            Turb(TurbNum).RotorNormVec = [cos(Turb(TurbNum).YawAngle); sin(Turb(TurbNum).YawAngle)];
            Turb(TurbNum).ThrustForceVec = [0; 0];
            Output_RotorAcc(TurbNum) = 0;
        end

        %% 更新速度历史（用于辐射力）
        if Time > 0
            new_buffer_ptr = mod(new_buffer_ptr, memory_length) + 1;
            new_vel_buffer(:, new_buffer_ptr) = Turb(TurbNum).VelVec;
        end

        % ============================================================
        % 水动力：直接使用 MATLAB Function1 输出的 WaveExcitationForce
        % 不在 StateDerivVec 内再次调用 test_cal_FK_diff。
        % 这样可避免同一时间步重复计算水动力、重复维护 persistent 状态，
        % 并保证 Simulink 中唯一的水动力源就是 MATLAB Function1。
        % ============================================================
        if size(WaveExcitationForce,1) >= 2 && size(WaveExcitationForce,2) >= TurbNum
            Turb(TurbNum).HydroForceVec = WaveExcitationForce(1:2,TurbNum);
        else
            Turb(TurbNum).HydroForceVec = [0; 0];
        end

        %% 连续时间水动力诊断：每 1 s 打印一次
            if Time >= 0 && (Time - last_hydro_diag_time >= 1.0 - 1e-9)
            fprintf('[HydroState] T=%8.3f: HydroForceVec=[% .6e,% .6e] N\n', ...
                Time, Turb(TurbNum).HydroForceVec(1), Turb(TurbNum).HydroForceVec(2));
            last_hydro_diag_time = Time;
        end

        %% Mooring line restoring force
        FLPosRelAnch1 = Turb(TurbNum).PosVec + Turb(TurbNum).FL1PosRelG - Turb(TurbNum).Anch1PosVec;
        FLPosRelAnch2 = Turb(TurbNum).PosVec + Turb(TurbNum).FL2PosRelG - Turb(TurbNum).Anch2PosVec;
        FLPosRelAnch3 = Turb(TurbNum).PosVec + Turb(TurbNum).FL3PosRelG - Turb(TurbNum).Anch3PosVec;
        Line1HorDist = norm(FLPosRelAnch1);
        Line2HorDist = norm(FLPosRelAnch2);
        Line3HorDist = norm(FLPosRelAnch3);
        if length(HorDistData) >= 2 && length(HorForceData) >= 2
            Line1HorForce = interp1(HorDistData, HorForceData, Line1HorDist, 'linear', 0);
            Line2HorForce = interp1(HorDistData, HorForceData, Line2HorDist, 'linear', 0);
            Line3HorForce = interp1(HorDistData, HorForceData, Line3HorDist, 'linear', 0);
        else
            Line1HorForce = 0;
            Line2HorForce = 0;
            Line3HorForce = 0;
        end
        if Line1HorDist > 0
            MoorForce1 = -Line1HorForce * FLPosRelAnch1 / Line1HorDist;
        else
            MoorForce1 = [0; 0];
        end
        if Line2HorDist > 0
            MoorForce2 = -Line2HorForce * FLPosRelAnch2 / Line2HorDist;
        else
            MoorForce2 = [0; 0];
        end
        if Line3HorDist > 0
            MoorForce3 = -Line3HorForce * FLPosRelAnch3 / Line3HorDist;
        else
            MoorForce3 = [0; 0];
        end
        Turb(TurbNum).MoorForceVec = MoorForce1 + MoorForce2 + MoorForce3;

        %% 二次粘性阻尼 (文献 Table XII)
        %B_surge = 1.25e6;   % N*s^2/m^2  (Surge)
        %B_sway = 0.95e6;    % N*s^2/m^2  (Sway)

        %ViscousDampingForce = -[B_surge * abs(Turb(TurbNum).VelVec(1)) * Turb(TurbNum).VelVec(1);
        %                         B_sway * abs(Turb(TurbNum).VelVec(2)) * Turb(TurbNum).VelVec(2)];

        LinDampCoeff = 1000000;  % N/(m/s)
        LinDampForceVec = -LinDampCoeff * Turb(TurbNum).VelVec;

        %TotalDampingForce = ViscousDampingForce + LinDampForceVec;
        TotalDampingForce =LinDampForceVec;

        %% Turbine acceleration
        Turb(TurbNum).ForceVec = Turb(TurbNum).ThrustForceVec + Turb(TurbNum).MoorForceVec + Turb(TurbNum).HydroForceVec + TotalDampingForce;

        Turb(TurbNum).AccelVec = Solver.TurbMotion * Turb(TurbNum).ForceVec / (Turb(TurbNum).Mass + Turb(TurbNum).AddedMass);

        %% Turbine power calculation
        if enableWind
            if modelselect == 1
                if preshut(TurbNum) == 0
                    Turb(TurbNum).Power = Turb(TurbNum).RotorSpeed * Turb(TurbNum).GeneTorque;
                    TurbPowerVec(TurbNum) = Turb(TurbNum).Power;
                else
                    TurbPowerVec(TurbNum) = 0;
                end
            else
                Turb(TurbNum).Power = 0.5 * Turb(TurbNum).PowerCoeff * Global.AirDensity * Turb(TurbNum).RotorArea * (norm(Turb(TurbNum).RelUpsWindVelVec)^3);
                TurbPowerVec(TurbNum) = Turb(TurbNum).Power;
            end
        else
            Turb(TurbNum).Power = 0;
            TurbPowerVec(TurbNum) = 0;
        end

        %% Calculate wake state derivatives
        if enableWind
            if strcmp(Solver.RotorMomModel, 'Jimenez')
                InitWakeVelMag = norm(Turb(TurbNum).RelUpsWindVelVec) * sqrt(abs(1 - Turb(TurbNum).ThrustCoeff));
                InitWakeAngle = -((cos(Turb(TurbNum).RelYawAngle))^2) * sin(Turb(TurbNum).RelYawAngle) * Turb(TurbNum).ThrustCoeff / 2;
                Turb(TurbNum).InitWakeVelVec = InitWakeVelMag * [
                    cos(InitWakeAngle + Turb(TurbNum).RelWindAngle);
                    sin(InitWakeAngle + Turb(TurbNum).RelWindAngle)];
            elseif strcmp(Solver.RotorMomModel, 'Bastankhah')
                InitWakeVelMag = norm(Turb(TurbNum).RelUpsWindVelVec) * sqrt(1 - Turb(TurbNum).ThrustCoeff);
                InitWakeAngle = -0.3 * Turb(TurbNum).RelYawAngle / cos(Turb(TurbNum).RelYawAngle) * (1 - sqrt(1 - Turb(TurbNum).ThrustCoeff * cos(Turb(TurbNum).RelYawAngle)));
                Turb(TurbNum).InitWakeVelVec = InitWakeVelMag * [
                    cos(InitWakeAngle + Turb(TurbNum).RelWindAngle);
                    sin(InitWakeAngle + Turb(TurbNum).RelWindAngle)];
            else
                error('Selected momentum model does not exist...');
            end

            for PointNum = 1:Turb(TurbNum).NumWakePoints
                PropVel = [1 0] * Global.FreeStreamWindVel - Turb(TurbNum).VelVec(1);
                WakeRecovGrad = (1/(pi/4*Turb(TurbNum).WakePointDia(PointNum)^2)) * (pi/2*Turb(TurbNum).WakePointDia(PointNum)*Global.WakeExpRate) * (Global.FreeStreamWindVel - (Turb(TurbNum).VelVec + Turb(TurbNum).WakePointVelVec(:, PointNum)));
                if PointNum == 1
                    WakeLocalVel(PointNum) = Turb(TurbNum).WakePointVelVec(2, PointNum) - PropVel * (Turb(TurbNum).WakePointPos_Y(PointNum + Solver.CentralDiff)) / Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff);
                    WakeLocalAccel(:, PointNum) = Global.FreeStreamWindAccel - Turb(TurbNum).AccelVec - PropVel * (Turb(TurbNum).WakePointVelVec(:, PointNum + Solver.CentralDiff) - Turb(TurbNum).InitWakeVelVec) / Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff) + WakeRecovGrad;
                    WakeLocalDiaROT(PointNum) = Global.WakeExpRate - PropVel * (Turb(TurbNum).WakePointDia(PointNum + Solver.CentralDiff) - Turb(TurbNum).RotorDia) / Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff);
                elseif PointNum == Turb(TurbNum).NumWakePoints
                    WakeLocalVel(PointNum) = Turb(TurbNum).WakePointVelVec(2, PointNum) - PropVel * (Turb(TurbNum).WakePointPos_Y(PointNum) - Turb(TurbNum).WakePointPos_Y(PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum) - Turb(TurbNum).WakePointPos_X(PointNum - 1));
                    WakeLocalAccel(:, PointNum) = Global.FreeStreamWindAccel - Turb(TurbNum).AccelVec - PropVel * (Turb(TurbNum).WakePointVelVec(:, PointNum) - Turb(TurbNum).WakePointVelVec(:, PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum) - Turb(TurbNum).WakePointPos_X(PointNum - 1)) + WakeRecovGrad;
                    WakeLocalDiaROT(PointNum) = Global.WakeExpRate - PropVel * (Turb(TurbNum).WakePointDia(PointNum) - Turb(TurbNum).WakePointDia(PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum) - Turb(TurbNum).WakePointPos_X(PointNum - 1));
                else
                    WakeLocalVel(PointNum) = Turb(TurbNum).WakePointVelVec(2, PointNum) - PropVel * (Turb(TurbNum).WakePointPos_Y(PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointPos_Y(PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointPos_X(PointNum - 1));
                    WakeLocalAccel(:, PointNum) = Global.FreeStreamWindAccel - Turb(TurbNum).AccelVec - PropVel * (Turb(TurbNum).WakePointVelVec(:, PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointVelVec(:, PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointPos_X(PointNum - 1)) + WakeRecovGrad;
                    WakeLocalDiaROT(PointNum) = Global.WakeExpRate - PropVel * (Turb(TurbNum).WakePointDia(PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointDia(PointNum - 1)) / (Turb(TurbNum).WakePointPos_X(PointNum + Solver.CentralDiff) - Turb(TurbNum).WakePointPos_X(PointNum - 1));
                end
            end

            WakeStateDerivMat(:, BaseNum + 1:BaseNum + Turb(TurbNum).NumWakePoints) = [
                WakeLocalVel(1:Turb(TurbNum).NumWakePoints);
                WakeLocalAccel(:, 1:Turb(TurbNum).NumWakePoints);
                WakeLocalDiaROT(1:Turb(TurbNum).NumWakePoints)];
        else
            Turb(TurbNum).InitWakeVelVec = [0; 0];
            WakeStateDerivMat(:, BaseNum + 1:BaseNum + Turb(TurbNum).NumWakePoints) = 0;
        end

        TurbPosStateDerivMat(:, TurbNum) = Turb(TurbNum).VelVec;
        TurbVelStateDerivMat(:, TurbNum) = Turb(TurbNum).AccelVec;

        BaseNum = BaseNum + Turb(TurbNum).NumWakePoints;
    end

    %% Build final state derivative vector
    StateDerivs = [
        reshape(TurbPosStateDerivMat, [], 1);
        reshape(TurbVelStateDerivMat, [], 1);
        reshape(WakeStateDerivMat, [], 1)];

    for i = 1:length(StateDerivs)
        if isnan(StateDerivs(i)) || isinf(StateDerivs(i))
            StateDerivs(i) = 0;
        end
    end

    %% Outputs
    for TurbNum = 1:Farm.NumTurb
        Output_PosX(TurbNum) = Turb(TurbNum).PosVec(1);
        Output_PosY(TurbNum) = Turb(TurbNum).PosVec(2);
        Output_Power(TurbNum) = Turb(TurbNum).Power;
        Output_ThrustC(TurbNum) = Turb(TurbNum).ThrustCoeff;
        Output_PowerC(TurbNum) = Turb(TurbNum).PowerCoeff;
        Output_RotorSpeedv(TurbNum) = Turb(TurbNum).RotorSpeed;
        windx(TurbNum) = Turb(TurbNum).UpsWindVelVec(1);
        if enableWind
            if Turb(TurbNum).UpsWindVelVec(1) >= 3
                shutdown(TurbNum) = 0;
            else
                shutdown(TurbNum) = 1;
                Output_Power(TurbNum) = 0;
            end
        else
            shutdown(TurbNum) = 0;
        end
    end

end