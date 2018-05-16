%% CIET Fault Detection Generic Script

% Christopher Poresky
% This script should be used to formulate a generic fault detection
% algorithm to be employed in ARCO-CIET fault detection experiments and
% simulations. It consists of various other sub-projects that come from CE
% 295 homework assignments, code developed for frequency response testing,
% and new code made for this application.

% Initiated with actual code 5/16/18

%% General setup
clear; close all;
fs = 15;    % Font size for plots


%% Parameter Identification

% Prepare data
 % File with experimental data
 filename = '2017-10-26_Frequency_Response_Heater.csv';
    % This is an example. The user will need to input the filename.
    % IMPORTANT: Need to add necessary info about file formatting

 % Desired data series (Boolean) [Power Output, Heater Inlet, Heater
 % Outlet, CTAH Outlet, CTAH Inlet, Heater Surface]
 desired_params = [1 1 1 0 0 1];
    % This is an example. The user will need to choose the parameters of
    % interest for the given situation.
    % IMPORTANT: Currently, this is made to work with DataExtract.m, which
    % is designed for CIET data of a certain format that may have changed
    % recently and may not remain the same as we move forward with
    % ARCO/CIET work.
    
 % Store time vector, data vectors, and name vectors
 [time,data_array,name_array] = DataExtract(filename,desired_params);
    % This is generic but requires the use of DataExtract.m, which is
    % currently NOT generic for all potentially-relevant CIET parameters.
    
 % Define data bounds for training
 data_start = 7566; % seconds
 data_stop = 8828; % seconds
    % This is an example. The user will need to define the start and end
    % times of the test.
    % IMPORTANT: Currently, the start and end times are defined by visual
    % inspection and recorded procedures. In the future, it may be better
    % to use a system that automates start and end time selection or uses
    % some more explicity criteria.
    
 [data_array_trim,time_trim] = DataTrim(data_array,time,...
     data_start,data_stop);
    % This is generic but it would be good to thoroughly inspect DataTrim.m
    % to confirm generality.
    
 % Create data matrix
 write_matrix = time_trim;
 for i = 1:length(data_array_trim)
     write_matrix = [write_matrix data_array_trim{i}];
 end
    % This is obviously not optimized computationally. It would be best to
    % improve efficiency sooner rather than later. The point here is to
    % create an array that has all data series and time series together.
    % Perhaps this step could be avoided if the time values and data values
    % were never separated in the first place.
 
 csvwrite('Param_ID_Data.csv',write_matrix);
    % This filename could change and use some kind of automated naming
    % convention so that multiple files can be generated successively
    % without needing the user to externally modify filenames, relocate, or
    % even delete things.

 % Load data and assign to matrix
 data = csvread('Param_ID_Data.csv'); 
 t = data(:,1);     % t    : time vector [sec]
 T_in = data(:,3);  % T_in : heater inlet temperature [deg C]
 T_F = data(:,4);   % T_F  : heater outlet temperature [deg C]
 T_S = data(:,5);   % T_S  : heater surface temperature [deg C]
 P = data(:,2);     % P    : heater power [W]
    % It's important that the argument for csvread.m is the same as the
    % argument for csvwrite.m above. As mentioned above, it could be
    % standardized to depend on the specific run. Additionally, all of the
    % variable assignments here are examples and case-dependent.

% Check persistence of excitation

 % Define parameter vectors
 phi = zeros(2,length(t),2);
    % This part depends on how many parameters you have. So, 2 here is an
    % example and the user will need to modify accordingly.
 phi(:,:,1) = [T_in - T_F, T_S - T_F]';
 phi(:,:,2) = [T_F - T_S, P]';
    % This is an example. The definition of the parameter vectors is a
    % primary user-input situation.
 
 t_end = t(end) - t(1);
    % This is generic.
    
 PE_mat = zeros(size(phi,1),size(phi,1),size(phi,3));
    % This is generic.
    
 phi_sq = zeros(size(phi,1),size(phi,1),size(phi,3),length(t));
 for i = 1:size(phi_sq,3)
     for k = 1:length(t)
         phi_sq(:,:,i,k) = phi(:,k,i) * phi(:,k,i)';
     end
     for m = 1:size(phi_sq,1)
         for n = 1:size(phi_sq,2)
             PE_mat(m,n,i) = 1/t_end * trapz(t,phi_sq(m,n,i,:));
         end
     end
 end
 
 PE_lam_min - zeros(1,size(phi,3));
 for i = 1:length(PE_lam_min)
     PE_lam_min(i) = min(min(eig(PE_mat(:,:,i))));
 end
 
 fprintf(1,'PE Level for Separated Parameters, T_bar_F : %1.4f\n',PE_lam_min(1));
 fprintf(1,'PE Level for Separated Parameters, T_bar_S : %1.4f\n',PE_lam_min(2));

% Plot data vs. time

% Feed data through gradient algorithm to identify parameters

% Examine parameter convergence times

% Build identified model

%% Validation

% Load validation data

% Simulate model and compare with validation data


%% Observer Design

% Define parameters

% Perform an observability analysis

% Load data and visualize

% Build open-loop observer (Luenberger)

% Find eigenvalues of the open-loop system (A matrix)

% Calculate the desired poles of the estimation error system (multiply
% eigenvalues of A matrix by some scalar between 1 and 10)

% Compute the observer gain using the place command

% Define the Luenberger state-space matrices

% Build the state-space model

% Define observer inputs, initial conditions

% Simulate observer response

% Build closed-loop observer (KF)