%Function to simulate rhmODE for WT hosts
% inputs: (1)Ki- total immune capacity, (2)immune- initial immune response, 
%         (3)phage-sensitive bacteria- initial bacteria conc., (4)antibiotic-sensitive bacteria- initial bacteria conc.
%         (5)phage- initial phage conc.
%         (6) antibiotic conc., - initial injected dose, (7) Antibiotic name

% outputs:(1)y- matrix of all of the population changes for bacteria,
%          phage and immune response.
%         (2)TB - vector of total bacteria population
%         (3)time - vector of time values over the simulation

function [y, TB, time ] = simRHM_WT_sensitivity(Ki, immune, Bp, Ba, phage, varargin)

    antibiotic = varargin{1}; % ug/ml, concentration of antibiotic
    p = varargin{2}; % structure containing param. values

    %initial conditions
    Bo = Bp;
    Ro = Ba;
    Po = 0;
    Io = immune;
    Ao = 0;
    tspan = 0:2;
    y0 = [Bo;Ro;Po;Io;Ao];

    % simulating diff eq until time to add phage (2hrs)
    options = odeset('Events',@myEventsFcn, 'RelTol', 1e-8);
    [t1,y1] = ode45(@rhmODE,tspan,y0,options,p);

    %----------------------------------------
    % Add Phage and Antibiotic dose time delay (2hrs)

    B = y1(end,1);
    R = y1(end,2);
    P = phage;
    I = y1(end,4);
    A = antibiotic;
    tspan2 = 2:96;
    yi = [B;R;P;I;A];
    p.pulse = antibiotic * p.theta;
    % simulating diff eq after phage and antibiotic addition
    [t2,y2] = ode45(@rhmODE,tspan2,yi,options,p);
    

    B = y2(end,1);
    R = y2(end,2);
    % Check whether after the addition of the first
    % dose of antibiotic or after several hours of simulation run
    % the Bacteria reach extinction.
    if B <= 1
        B = 0;
    end
    if R <= 1
        R = 0;
    end
    P = y2(end,3);
    I = y2(end,4);
    
    %----------------------------------------
    % continue simulation after susceptibles or resistants die

    check = 0;
    currentTime = t2(end);
    if currentTime < 96 % Bacterial pop died before end of simulation
        check = 1;
        B2 = y2(end,1);
        if B2 <= 1
            B2 = 0;
        end
        R2 = y2(end,2);
        if R2 <= 1
            R2 = 0;
        end
        P2 = y2(end,3);
        I2 = y2(end,4);
        A2 = y2(end,5);
        tspan3 = currentTime:0.05:96;
        yii = [B2;R2;P2;I2;A2];

        % simulating diff eq
        [t3,y3] = ode45(@rhmODE,tspan3,yii,options,p);
    end
    %----------------------------------------
    % continue simulation after susceptibles or resistants die

    % Check if run completed without bacteria dying
    if check == 1
        currentTime2 = t3(end);
    else
        currentTime2 = 96;
    end
    if currentTime2 <= 95
        B3 = y3(end,1);
        if B3 <= 1
            B3 = 0;
        end
        R3 = y3(end,2);
        if R3 <= 1
            R3 = 0;
        end
        P3 = y3(end,3);
        I3 = y3(end,4);
        A3 = y3(end, 5);
        tspan4 = currentTime2:96;
        yiii = [B3;R3;P3;I3;A3];

        % simulating diff eq
        [t4,y4] = ode45(@rhmODE,tspan4,yiii,options,p);
        
        currentTime3 = t4(end);
        
        if currentTime3 <= 95
               B4 = y4(end,1);
            if B4 <= 1
                B4 = 0;
            end
            R4 = y4(end,2);
            if R4 <= 1
                R4 = 0;
            end
            P4 = y4(end,3);
            I4 = y4(end,4);
            A4 = y4(end,5);
            tspan5 = currentTime3:96;
            yiiii = [B4;R4;P4;I4;A4];

            % simulating diff eq
            [t5, y5] = ode45(@rhmODE,tspan5,yiiii,options,p);
            t4 = [t4; t5];
            y4 = [y4; y5];                    
        end

        time = [t1(1:end-1); t2; t3; t4];
        y = [y1(1:end-1, :); y2; y3; y4];
    elseif check == 1
        time = [t1(1:end-1); t2; t3];
        y = [y1(1:end-1,:); y2; y3];
    else
        time = [t1(1:end-1); t2];
        y = [y1(1:end-1,:); y2];
    end
    S = y(:,1);
    R = y(:,2);
    TB = S + R; 
end

