function mainFunction
    clear all;
    clc;

    global recordedData; % Access the global variable
    
    % Run the EM GUI
    electromigration_GUI;
    
    % Access the recorded data after the GUI closes
    disp('Recorded Data:');
    disp(recordedData);
    
    %% Simulation Parameters

    % Spatial parameters based on geometry
    A=extract(recordedData);
    L_list = A(1,:)*1e-6;                        % lengths of segments
    W_list = A(2,:)*1e-6;                        % widths of segments
    H_list = A(3,:)*1e-6;                        % heights of segments
    L_divs = 10;
    j_list = A(4,:)*1e9;                         % current densities of (m + 1 + n) segments

    % Simulation time parameters
    T=86400;
    t_step=8640;

    % Partial infinite sum parameters
    N = 3;       % no. of terms in partial sum
    start_n = 1; % n=1 start point of partial sum

    % Temperature and Compact Model parameters
    temp=recordedData.T;
    
    % Compact Model parameters
    A = model_params(recordedData.metal,recordedData.confinement, temp);
    alpha0 = A(1);
    beta = A(2);
    sigma_lim = A(3);
    gamma = A(4);
    
    alpha = alpha0;

    %% INTERCONNECT NETWORK TYPE--1
    if recordedData.f==1
        h=msgbox('...','Progress...');
        set(h,'Position', [300, 300, 500, 300]);
        
        % Geometrical parameters
        m=recordedData.m;                             % number of segments on left of central segment
        n=recordedData.n;                             % number of segements on right of central segment
        
        % Current injection
        I_inj = recordedData.I;                      % injection current magnitude
        
        I_inj = I_inj * (randi([0, 1]) * 2 - 1);     % injection current directionality (injection / withdrawl)
        inj_loc = randi([0, 1]);                     % injection location
        
        if inj_loc==0      % left-central node injection
            inj_L = I_inj;
            inj_R = 0;
        elseif inj_loc==1  % right-central node injection
            inj_L = 0;
            inj_R = I_inj;
        end
        
        % Updating current densities
        A_list=W_list.*H_list;
        I_list=j_list.*A_list;
        
        if inj_loc==0      % left-central node injection
            V_L=(inj_L +sum(I_list(1:m+1)))/sum(A_list(1:m+1)./L_list(1:m+1));
            j_list(1:m+1)=V_L./L_list(1:m+1);
        elseif inj_loc==1  % right-central node injection
            V_R=(inj_R +sum(I_list(m+1:m+n+1)))/sum(A_list(m+1:m+n+1)./L_list(m+1:m+n+1));
            j_list(m+1:m+n+1)=V_R./L_list(m+1:m+n+1);
        end
        
        
        tic; % start timer
        
        %% Stress parameters calculation
        set(findobj(h, 'Type', 'text'), 'String', 'Calculating stress parameters...\n', 'FontSize', 14);
        
        % c_i,1 calculation
        
        c1_m_list=beta*j_list(1:m);                % left terminal nodes
        c1_n_list=beta*j_list(m+2:length(j_list)); % right terminal nodes
        c1_central=(1/(W_list(m+1)*H_list(m+1)))*(sum(beta*W_list(1:m+1).*H_list(1:m+1).*j_list(1:m+1))-sum(W_list(1:m).*H_list(1:m).*c1_m_list));  % left-central node
        c1_central=(1/(W_list(m+1)*H_list(m+1)))*(sum(beta*W_list(m+1:length(W_list)).*H_list(m+1:length(H_list)).*j_list(m+1:length(j_list)))-sum(W_list(m+2:length(W_list)).*H_list(m+2:length(W_list)).*c1_n_list));  % right-central node
        
        % c_i,2 calculation
        
        syms c2_ [1 m+n+1]
        syms X
        
        equations=[c2_(m)==X];
        
        for i=1:m
            equations=[equations c1_m_list(i)*L_list(i)+c2_(i) == c2_(m+1)];
        end
        
        for i=1:n
            equations = [equations c2_(m+1+i) == c1_central*L_list(m+1)+c2_(m+1)];
        end
        
        %pretty(equations)
        c2_sol=struct2array(vpasolve(equations,c2_));
        phi_list=c2_sol-X;
        
        sum_left = sum( (W_list(1:m).*H_list(1:m)).*(0.5.*c1_m_list.*(L_list(1:m).^2) + phi_list(1:m).*L_list(1:m)) );
        sum_center = (W_list(m+1)*H_list(m+1)).*(0.5.*c1_central*(L_list(m+1)^2) + phi_list(m+1).*L_list(m+1));
        sum_right = sum( (W_list(m+2:length(W_list)).*H_list(m+2:length(H_list))).*(0.5.*c1_n_list.*(L_list(m+2:length(L_list)).^2) + phi_list(m+2:length(phi_list)).*L_list(m+2:length(L_list))) );
        
        X_sol=-(sum_left+sum_center+sum_right)/sum(W_list.*H_list.*L_list);
        
        c2_sol=phi_list+X_sol;
        
        t_calc=toc % Time taken for parameter calculation
        
        set(findobj(h, 'Type', 'text'), 'String', 'EM stress parameters calculated.\n', 'FontSize', 14);
        
        % %% Stress Evaluation
        % set(findobj(h, 'Type', 'text'), 'String', sprintf("Evaluating EM stress at t = %d seconds timestep...\n",T), 'FontSize', 14);
        % 
        % m_seg_stress_list = [];
        % for i = 1:m
        %     m_seg_stress_list = [m_seg_stress_list stresscalc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, T)];
        % end
        % 
        % central_seg_stress = stresscalc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), L_divs, alpha, beta, start_n, N, T);
        % 
        % n_seg_stress_list = [];
        % for i = 1:n
        %     n_seg_stress_list = [n_seg_stress_list stresscalc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), L_divs, alpha, beta, start_n, N, T)];
        % end
        % 
        % t_EM_eval=toc-t_calc % EM Stress Evaluation Time
        % 
        % set(findobj(h, 'Type', 'text'), 'String', sprintf("EM stress evaluated at t = %d seconds timestep.\n",T), 'FontSize', 14);
        
        %% Voiding Location and Voiding Time Calculation
        set(findobj(h, 'Type', 'text'), 'String', 'Calculating voiding locations and voiding timesteps...\n', 'FontSize', 14);
        
        m_void_params_list=[];
        for i =1:m
            m_void_params_list=[m_void_params_list void_calc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
        end
        
        central_void_params = void_calc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), alpha, beta, start_n, N, sigma_lim,T,t_step);
        
        n_void_params_list=[];
        for i =1:n
            n_void_params_list = [n_void_params_list void_calc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
        end
        
        t_total=toc; % Total Time taken
        t_void_calc=t_total-t_calc %-t_EM_eval % Time taken for Voiding Location & Voiding Time Calculation
        t_total
        
        set(findobj(h, 'Type', 'text'), 'String', 'Voiding Locations and Voiding Timesteps calculated.\n', 'FontSize', 14);
        
        
        fprintf("Storing data for this iteration...\n");
        %data(1,:)={double(m),double(n),mat2str(L_list),mat2str(W_list),mat2str(H_list),mat2str(j_list),mat2str([c1_m_list c1_central c1_n_list]),mat2str(double(c2_sol)),inj_loc,I_inj,double(t_calc),mat2str(double(m_seg_stress_list)), mat2str(double(central_seg_stress)), mat2str(double(n_seg_stress_list)),double(t_EM_eval),mat2str(double(m_void_params_list)), mat2str(double(central_void_params)), mat2str(double(n_void_params_list)),double(t_void_calc),double(t_total)};
        fprintf("Data stored for this iteration.\n");

        close(h)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Plotting
        % Calculate the number of rows in the grid (maximum of m or n)
        rows = max(m, n);
        
        % Create figure
        figure;
        
        
        % Define the height of each subplot
        subplot_height_m = 0.8 / m; % reduce total height to leave space below
        subplot_height_n = 0.8 / n; % reduce total height to leave space below
        
        % Plot subplots in the first column
        for i = 1:m
            pos = [0.1, 0.1 + (m-i) * subplot_height_m, 0.25, subplot_height_m - 0.12]; % [left, bottom, width, height]
            subplot('Position', pos);
            hold on
            %f=0;
            for t=0:t_step:T
                x=linspace(0,L_list(i),L_divs);
                y=stresscalc(c1_m_list(i),c2_sol(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, t);
                plot(x,y);
                idx_above = find(y >= sigma_lim);
                idx_below = find(y <= -sigma_lim);
                plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
                plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
            end
            hold off
            title(sprintf('Left Segment %d\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', i, L_list(i)/1e-6,W_list(i)/1e-6,H_list(i)/1e-6,j_list(i)));
            xlim([0 L_list(i)])
        end
        
        % Plot the subplot in the second column (centered vertically)
        pos = [0.4, 0.3, 0.25, 0.4]; % [left, bottom, width, height]
        subplot('Position', pos);
        hold on
        %f=0;
        for t=0:t_step:T
            x=linspace(0,L_list(m+1),L_divs);
            y=stresscalc(c1_central,c2_sol(m+1),j_list(m+1), L_list(m+1), L_divs, alpha, beta, start_n, N, t);
            plot(x,y);
            idx_above = find(y >= sigma_lim);
            idx_below = find(y <= -sigma_lim);
            plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
            plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
        end
        hold off
        title(sprintf('Central Segment\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', L_list(m+1)/1e-6,W_list(m+1)/1e-6,H_list(m+1)/1e-6,j_list(m+1)));
        xlim([0 L_list(m+1)])
        
        % Plot subplots in the third column
        for i = 1:n
            pos = [0.7, 0.1 + (n-i) * subplot_height_n, 0.25, subplot_height_n - 0.12]; % [left, bottom, width, height]
            subplot('Position', pos);
            hold on
            %f=0;
            for t=0:t_step:T
                x=linspace(0,L_list(m+1+i),L_divs);
                y=stresscalc(c1_n_list(i),c2_sol(m+1+i),j_list(m+1+i), L_list(m+1+i), L_divs, alpha, beta, start_n, N, t);
                plot(x,y);
                idx_above = find(y >= sigma_lim);
                idx_below = find(y <= -sigma_lim);
                plot(x(idx_above), y(idx_above), 'ro', 'MarkerFaceColor', 'r');
                plot(x(idx_below), y(idx_below), 'ro', 'MarkerFaceColor', 'r');
            end
            hold off
            title(sprintf('Right Segment %d\nLength: %d um, Width: %.4f um, Height: %.4f um\nj = %d A/m^2', i, L_list(m+1+i)/1e-6,W_list(m+1+i)/1e-6,H_list(m+1+i)/1e-6,j_list(m+1+i)));
            xlim([0 L_list(m+1+i)])
        end


    %% INTERCONNECT NETWORK TYPE--2    
    else
        h=msgbox('...','Progress...');
        set(h,'Position', [300, 300, 500, 300]);

        % Geometrical Parameters
        n = recordedData.n;                            % height of tree
        
        % Current injection
        I_inj = recordedData.I;                        % injection current magnitude
        
        I_inj = I_inj * (randi([0, 1]) * 2 - 1);       % injection current directionality (injection / withdrawl)
        inj_loc = randi([0, 1]);                       % injection location
        
        tic; % start timer
        
        %% Stress parameters calculation
        set(findobj(h, 'Type', 'text'), 'String', 'Calculating stress parameters...\n', 'FontSize', 14);       
        
        % c_i,1 calculation
        
        c1_start=beta*j_list(1);                       % left terminal node
        c1_end_list=beta*j_list(end-2^(n-1)+1:end);    % right terminal nodes
        
        syms c1_ [1 2^n-1]
        
        c1_(1)=c1_start;
        c1_(end-2^(n-1)+1:end)=c1_end_list;
        
        equations=[];
        
        % Network-segment mapping
        A = adjacency_matrix_binTree(n);
        seg_map=generate_segment_map(2^n-1);
        
        
        for i=2:2^(n-1)
            segments_list=connected_segments(A,i);
            L=[];
            for j=1:length(segments_list)
                segment=segments_list(j,:);
                L=[L get_segment_number(segment(1),segment(2),seg_map)];
            end
            sum_l=0;
            sum_r=0;
            for k=1:length(L)
                sum_l = sum_l + W_list(L(k))*H_list(L(k))*c1_(L(k));
                sum_r = sum_r + W_list(L(k))*H_list(L(k))*j_list(L(k));
            end
            equations=[equations sum_l==beta*sum_r];
        end
        c1_(2:2^(n-1)-1)=struct2array(vpasolve(equations(1:end-1),c1_(2:2^(n-1)-1)));
        
        c1_=double(c1_);
        
        % c_i,2 calculation
        
        syms c2_ [1 2^n-1]
        syms X
        
        equations=[c2_(1)==X];
        
        for i=2:2^(n-1)
            segments_list=connected_segments(A,i);
            L=[];
            for j=1:length(segments_list)
                segment=segments_list(j,:);
                L=[L get_segment_number(segment(1),segment(2),seg_map)];
            end
            equations=[equations c1_(L(1))*L_list(L(1))+c2_(L(1))==c2_(L(2)) c1_(L(1))*L_list(L(1))+c2_(L(1))==c2_(L(3))];
        end
        
        c2_=struct2array(vpasolve(equations,c2_));
        phi_list=c2_-X;
        
        X_sol=-sum(W_list.*H_list.*(0.5.*c1_.*(L_list.^2) + phi_list.*L_list))/sum(W_list.*H_list.*L_list);
            
        c2_=phi_list+X_sol;
        c2_=double(c2_);
        
        t_calc=toc % Time taken for parameter calculation
        
        set(findobj(h, 'Type', 'text'), 'String', 'EM stress parameters calculated.\n', 'FontSize', 14);       
       
        % %% Stress Evaluation
        % set(findobj(h, 'Type', 'text'), 'String', sprintf ("Evaluating EM stress at t = %d seconds timestep...\n",T), 'FontSize', 14);
        % 
        % stress_list = [];
        % for i = 1:2^n-1
        %     stress_list = [stress_list stresscalc(c1_(i),c2_(i),j_list(i), L_list(i), L_divs, alpha, beta, start_n, N, T)];
        % end
        % stress_list;
        % 
        % t_EM_eval=toc-t_calc; % EM Stress Evaluation Time
        % 
        % set(findobj(h, 'Type', 'text'), 'String', sprintf("EM stress evaluated at t = %d seconds timestep.\n",T), 'FontSize', 14);

        %% Voiding Location and Voiding Time Calculation
        set(findobj(h, 'Type', 'text'), 'String', 'Calculating voiding locations and voiding timesteps...\n', 'FontSize', 14);
        
        void_params_list=[];
        for i =1:2^n-1
            void_params_list=[void_params_list void_calc(c1_(i),c2_(i),j_list(i), L_list(i), alpha, beta, start_n, N, sigma_lim,T,t_step)];
        end
        
        void_params_list;
        
        t_total=toc; % Total Time taken
        t_void_calc=t_total-t_calc %-t_EM_eval; % Time taken for Voiding Location & Voiding Time Calculation
        t_total
        
        set(findobj(h, 'Type', 'text'), 'String', 'Voiding Locations and Voiding Timesteps calculated.\n', 'FontSize', 14);
        
    
        fprintf("Storing data for this iteration...\n");
        %data(sim_n+1,:)={double(n),double(2^n - 1),mat2str(L_list),mat2str(W_list),mat2str(H_list),mat2str(j_list),mat2str([c1_]),mat2str(c2_),double(t_calc),mat2str(double(stress_list)),double(t_EM_eval),mat2str(double(void_params_list)),double(t_void_calc),double(t_total)};
        fprintf("Data stored for this iteration.\n");

        close(h)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        %% Plotting
        generate_plots(n, c1_, c2_, j_list, L_list, L_divs, alpha, beta, sigma_lim, start_n, N, t_step, T)

        
    end



% end of mainFunction
end


%% Functions
%% Compact Model Parameter Function
function A=model_params(metal,confinement,T)
    if strcmp(confinement,'SiO2')
        a_Conf=0.05e-5;
    elseif strcmp(confinement,'Si3N4')
        a_Conf=0.33e-5;
    end

    if metal=='Cu'
        alpha0 = ((24e-5*1.09e11*1.18e-29)/(T*1.38e-23))*exp(-3.25e-19/(T*1.38e-23));
        beta = ((1.602e-19)*(5))*(1.7e-8)/(1.18e-29); %(eZ*rho_Cu)/omega
        sigma_lim = 500e6-(1.09e11)*(1.65e-5-a_Conf)*(473-T);
        gamma = -(3.2e-21)/(T*1.38e-23);
    else
        alpha0 = ((12e-5*0.7e11*1.66e-29)/(T*1.38e-23))*exp(-2.35e-19/(T*1.38e-23));
        beta = ((1.602e-19)*(10))*(2.7e-8)/(1.66e-29); %(eZ*rho_Cu)/omega
        sigma_lim = 400e6-(0.7e11)*(2.31e-5-a_Conf)*(573-T);
        gamma = -(2.4e-21)/(T*1.38e-23);
    end
    A=[alpha0,beta,sigma_lim,gamma];
end
%% EM Stress Calculation Function
function result=stresscalc(c1,c2,j, L, L_divs, alpha, beta, start_n, N, t)
    x = linspace(0, L, L_divs); % entire interconnect check
    % FAST
    % x = [x(1),x(end)]; % interconnect endpoint checks only
    result=zeros(1,length(x));
    for i=1:length(x)
        syms n
        f1=(((-1)^n-1)/(n^2))*cos(n*pi*x(i)/L)*exp((-alpha*n^2*pi^2*t)/(L^2));
        f2=(((-1)^n-1)/(n^2))*sin(n*pi*x(i)/L)*exp((-alpha*n^2*pi^2*t)/(L^2));
        result(i)=c1*x(i) + c2 - ((2*c1*L)/pi^2)*symsum(f1, n, start_n, N) + (2/pi)*( c2 + (c1*L)/2 )*symsum(f2, n, start_n, N);
    end
end

%% Void Location & Void Time Calculation Function
function result=void_calc(c1, c2, j, L, alpha, beta, start_n, N, sigma_lim, T, t_step)
    y0 = c2;
    yL = c1 * L + c2;
    if abs(y0) > abs(yL)
        max_point_x = 0;
        max_point_y = y0;
    else
        max_point_x = L;
        max_point_y = yL;
    end

    loc=max_point_x;

    syms void_t
    syms n
    f1=(((-1)^n-1)/(n^2))*cos(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2));
    f2=(((-1)^n-1)/(n^2))*sin(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2));
    f = @(t) double(subs(c1*loc + c2 - ((2*c1*L)/pi^2)*symsum((((-1)^n-1)/(n^2))*cos(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2)), n, start_n, N) + (2/pi)*( c2 + (c1*L)/2 )*symsum((((-1)^n-1)/(n^2))*sin(n*pi*loc/L)*exp((-alpha*n^2*pi^2*void_t)/(L^2)), n, start_n, N),void_t,t));
    
    void_time=T;
    for t=0:t_step/10:T
        if f(t)>=sigma_lim
            void_time=t;
            break;
        end
    end
    result=[loc,void_time];
end


function A = adjacency_matrix_binTree(n)
 % Calculate the total number of nodes
    total_nodes = 2^n - 1;
    
    % Initialize the adjacency matrix
    A = zeros(total_nodes);
    
    % Iterate through each node and connect it to its children
    for i = 1:total_nodes
        left_child = 2*i;
        right_child = 2*i + 1;
        
        % If left child exists, create the connection
        if left_child <= total_nodes
            A(i, left_child) = 1;
            A(left_child, i) = 1;
        end
        
        % If right child exists, create the connection
        if right_child <= total_nodes
            A(i, right_child) = 1;
            A(right_child, i) = 1;
        end
    end
A1=zeros(size(A)+1);
A1(1,2)=1;A1(2,1)=1;
A1(2:end,2:end)=A;
A=A1;
end

function segments = connected_segments(A, i)
    % A: Adjacency matrix (nxn matrix)
    % i: The node index (integer)
    
    % Find all nodes connected to node i
    connected_nodes = find(A(i, :) == 1);  % Find indices where A(i,j) == 1
    
    % Create segments list as pairs of (i, j)
    segments = [repmat(i, length(connected_nodes), 1), connected_nodes'];
    for i=1:length(segments)
        if segments(i,1)>segments(i,2)
            L=segments(i,1);
            segments(i,1)=segments(i,2);
            segments(i,2)=L;
        end
    end
end

function segment_map = generate_segment_map(total_nodes)
    % total_nodes: the total number of nodes in the binary tree
    % segment_map: a matrix of pairs (parent, child) and their corresponding segment number

    segment_map = [];  % Initialize an empty matrix to store segment mappings
    segment_count = 1; % Start segment numbering from 1

    % Iterating over each node, starting from node 1
    for node = 1:total_nodes
        left_child = 2 * node;
        right_child = 2 * node + 1;

        % Checking if left child exists
        if left_child <= total_nodes
            segment_map = [segment_map; node, left_child, segment_count];  % Add the segment mapping
            segment_count = segment_count + 1;  % Increment segment number
        end

        % Check if right child exists
        if right_child <= total_nodes
            segment_map = [segment_map; node, right_child, segment_count];  % Add the segment mapping
            segment_count = segment_count + 1;  % Increment segment number
        end
    end
    segment_map=[[1 2 1];segment_map+1];
end

function k = get_segment_number(i, j, segment_map)
    % i: parent node
    % j: child node
    % segment_map: matrix containing all (i, j, k) mappings
    % k: the segment number corresponding to (i, j)

    % Finding the row in segment_map where the first two columns are i and j
    index = find(segment_map(:,1) == i & segment_map(:,2) == j);

    if ~isempty(index)
        k = segment_map(index, 3);  % Return the segment number k
    else
        k = -1;  % Return -1 if the segment (i,j) is not found
    end
end


function generate_plots(n, c1_, c2_, j_list, L_list, L_divs, alpha, beta, sigma_lim, start_n, N, t_step, T)
    % n: the number of layers
    f=0;
    figure_counter = 1;  % Counter for figure numbers

    for i = 1:n
        num_plots = 2^(i-1);  % Number of plots in the current layer
        figure(figure_counter);  % Create a new figure window
        figure_counter = figure_counter + 1;

        % Number of rows and columns for subplot layout
        rows = ceil(sqrt(num_plots));
        cols = ceil(num_plots / rows);

        for j = 1:num_plots
            subplot(rows, cols, j);  % Create a subplot
            hold on
            for t=0:t_step:T
                x=linspace(0,L_list(f+1),L_divs);
                y=stresscalc(c1_(i),c2_(f+1),j_list(f+1), L_list(f+1), L_divs, alpha, beta, start_n, N, t);
                plot(x,y);
            end
            yline(sigma_lim, '--r', ['y = ', num2str(sigma_lim)], 'LineWidth', 1.5);  % Add the y=constant line
            title(['Segment ', num2str(f+1)]);
            hold off;
            f=f+1;
        end

        % Adjust layout for better visibility
        sgtitle(['Layer ', num2str(i), ' with ', num2str(num_plots), ' segments']);
    end
end



%% GUI to spatial parameter extraction function
function A=extract(recData)
    f=recData.f;
    if f==1
        L_values = [str2double(recData.segmentValues(:,1))' str2double(recData.centralValues{1}) str2double(recData.additionalSegmentValues(:,1))'];
        W_values = [str2double(recData.segmentValues(:,2))' str2double(recData.centralValues{2}) str2double(recData.additionalSegmentValues(:,2))'];
        H_values = [str2double(recData.segmentValues(:,3))' str2double(recData.centralValues{3}) str2double(recData.additionalSegmentValues(:,3))'];
        j_values = [str2double(recData.segmentValues(:,4))' str2double(recData.centralValues{4}) str2double(recData.additionalSegmentValues(:,4))'];
    else
        L_values=[str2double(recData.segmentValues(:,1))'];
        W_values=[str2double(recData.segmentValues(:,2))'];
        H_values=[str2double(recData.segmentValues(:,3))'];
        j_values=[str2double(recData.segmentValues(:,4))'];
    end
    A=[L_values;W_values;H_values;j_values];
end