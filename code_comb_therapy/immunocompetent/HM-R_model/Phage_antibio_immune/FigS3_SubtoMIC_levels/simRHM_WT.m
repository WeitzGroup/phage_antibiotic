%Function to simulate rhmODE for Immunocompetent hosts
% inputs: (1)Ki- total immune capacity, (2)immune- initial immune response, 
%         (3)phage-sensitive bacteria- initial bacteria conc., (4)antibiotic-sensitive bacteria- initial bacteria conc.
%         (5)phage- initial phage conc.
%         (6) antibiotic conc., - initial injected dose, (7) Antibiotic name

% outputs:(1)y- matrix of all of the population changes for bacteria,
%          phage and immune response.
%         (2)TB - vector of total bacteria population
%         (3)time - vector of time values over the simulation

function [y, TB, time ] = simRHM_WT(Ki, immune, Bp, Ba, phage, varargin)

    if nargin == 5
        antibiotic = 0;
        anti_name = 'CP';
    else
        antibiotic = varargin{1};
        anti_name = varargin{2};
    end
    %---------parameters--------------
    % susceptible bacteria growth rate
    p.r = 0.75;
    % resistant bacteria growth rate
    p.rp = 0.675;
    % total bacteria carrying capacity
    p.Kc = 1e10;
    % nonlinear adsorption rate of phage:
    p.phi = 5.4e-8;
    % power law exponent in phage infection:
    p.g = 0.6;
    % immune response killing rate parameter:
    p.ep = 8.2e-8;
    % bacterial conc. at which immune response is half as effective:
    p.Kd = 4.1e7;
    % burst size of phage:
    p.beta = 100;
    % decay rate of phage:
    p.w = 0.07;
    % maximum growth rate of immune response:
    p.a = 0.97;
    % max capacity of immune response:
    p.Ki = Ki;
    % conc. of bacteria at which imm resp growth rate is half its maximum:
    p.Kn = 1e7;
    % probability of emergence of phage-resistant mutation per cell division
    p.m = 2.85e-8; %
    % probability of reversible mutation probability (Phage resistant -> Phage sensitive)
    p.m2 = 2.85e-8;
    
    if strcmp('CP', anti_name) % Ciprofloxacin parameters
        % Maximum antibiotic-mediated kill rate
        p.kkill = 18.5;
        % EC_50, antibiotic concentration when kkill/2
        p.ec = 0.3697;
        % Steepness of the sigmoid relationship between Antibiotic and death
        % rate
        p.H = 1;
        % Elimination rate of antibiotic from blood
        p.theta = 0.53;
        % Antibiotic pulse such that the Antibiotic conc. is at equilibrium
        % A = I - theta*A
        % A_star = I/theta
        % I = A_star*theta
        p.pulse = 0 * p.theta;
    
    else % Ceftazidime parameters
        % Maximum antibiotic-mediated kill rate
        p.kkill = 4.15;
        % EC_50, antibiotic concentration when kkill/2
        p.ec = 1.22;
        % Steepness of the sigmoid relationship between Antibiotic and death
        % rate
        p.H = 3.96;
        % Elimination rate of antibiotic from blood
        p.theta = 1.96;
        % Antibiotic pulse such that the Antibiotic conc. is at equilibrium
        p.pulse = 0 * p.theta;
    end
    

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
    tspan2 = 2:300;
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
    if currentTime < 300 % Bacterial pop died before end of simulation
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
        tspan3 = currentTime:300;
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
        currentTime2 = 300;
    end
    if currentTime2 <= 299
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
        tspan4 = currentTime2:300;
        yiii = [B3;R3;P3;I3;A3];

        % simulating diff eq
        [t4,y4] = ode45(@rhmODE,tspan4,yiii,options,p);
        
        currentTime3 = t4(end);

        if currentTime3 <= 299
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
            tspan5 = currentTime3:300;
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

