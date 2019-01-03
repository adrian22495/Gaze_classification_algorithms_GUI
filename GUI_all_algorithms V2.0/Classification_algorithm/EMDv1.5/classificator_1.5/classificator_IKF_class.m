% This code is based on Corey Holland's code.
% Department of Computer Science, Texas State University-San Marcos,
% San Marcos
% ch@txstate.edu

classdef classificator_IKF_class <  eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
                                    eye_records_class & ...                             % Basic class for placing eye tracker data
                                    eye_tracker_raw_data_converter_ETU_degree & ...     % Convertor between ETU and degress in data
                                    eye_tracker_raw_data_filter_class & ...             % Eye tracker data filtering by range of degrees
                                    classificator_merge_class & ...                     % Creates sequences of eye movements
                                    classificator_saccade_amplitude_filter_class & ...  % Filtered saccades based theire amplitude
                                    classificator_datafile_output_class & ...           % Output sequences to the files
                                    classificator_get_percentage_class & ...            % Calculate percentage of movements of every type
                                    classificator_enumerations_class & ...              % Basic enumerations definitions
                                    classificator_time_rate_class & ...                 % Time step and sample rate definition
                                    handle
    % This is implementation of I-KF classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup chi_threshold property calling set.chi_threshold method
    % 3. Setup window_size property calling set.window_size method
    % 4. Setup deviation property calling set.deviation method
    % 5. Call classify method to make classification
    
    properties (Hidden)
% Publicly accessible properties
        chi_threshold;      % Chi threshold for difference between saccade/threshold 
        window_size;        % Window size for chi-comparison
        deviation;          % Deviation value
    end
    
    methods
% Classification function
        function classify(obj)
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Begin data classification with I-KF classifier in :',datestr(now),'\n'));
            end
            obj.calculate_delta_t();
% Reset variables
            KF_K = 0;                          % This is Kalman filter's variables
            KF_x = [0; 0];                     % For any details look for Kalman filter information
            KF_y = [0; 0];                     % ...
            KF_P = [1 0; 0 1;];                % ...
            KF_A = [1 obj.delta_t_sec; 0 1];   % ...
            KF_H = [1 0];                      % ...
            KF_Q = [0.5 0; 0 0.5];             % ...
            KF_R = 0.5;                        % ...

            KF_predicted = zeros(length(obj.eye_records),1); % This array holds the predicted speeds
            KF_measured = zeros(length(obj.eye_records),1);  % This array holds the predicted speeds
            KF_x_velocity_degree = zeros(length(obj.eye_records),1);
            KF_y_velocity_degree = zeros(length(obj.eye_records),1);
% Calculate absolute degree velocity of our records
            KF_x_velocity_degree( 2:end ) =(    obj.eye_records( 2:end,obj.X_COORD ) - ...
                                                obj.eye_records( 1:end-1,obj.X_COORD) ) / obj.delta_t_sec;
            KF_y_velocity_degree( 2:end ) =(    obj.eye_records( 2:end,obj.Y_COORD ) - ...
                                                obj.eye_records( 1:end-1,obj.Y_COORD ) ) / obj.delta_t_sec;
% First point is a special case
            KF_x_velocity_degree(1) = 0;
            KF_y_velocity_degree(1) = 0;
            obj.eye_records(:,obj.VELOCITY) = sqrt( KF_x_velocity_degree.^2 + KF_y_velocity_degree.^2 );            
% First point is a special case
            obj.eye_records(1,obj.VELOCITY) = 0;

% Generating predicteded velocities using Kalman filter
            for i=1:length(obj.eye_records)
% Discarding the noise
                if ( obj.eye_records(i,obj.VALIDITY) == obj.DATA_INVALID || isnan(obj.eye_records(i,obj.VELOCITY)) )
                    obj.eye_records(i,obj.MOV_TYPE) = obj.NOISE_TYPE;
                else
% Prediction step of Kalman filter
                    KF_x = KF_A * KF_x;
                    KF_y = KF_A * KF_y;
                    KF_P = KF_A * KF_P * KF_A.' + KF_Q;
% Calculate predicteded velocity, storing current predicteded and measured velocities
                    KF_predicted(i) =   sqrt(KF_x(2).^2 + KF_y(2).^2);
                    KF_measured(i) =    obj.eye_records(i,obj.VELOCITY);
% Update step of Kalman filter
% Cory code
                    KF_K = (KF_P * KF_H.') / (KF_H * KF_P * KF_H.' + KF_R);
% Sam's code
%                    KF_K = KF_P.' * KF_H.' / (KF_H * KF_P.' * KF_H.' + KF_R);
                    KF_x = KF_x + KF_K * (KF_H * [obj.eye_records(i,obj.X_COORD); KF_x_velocity_degree(i)] - KF_H * KF_x);
                    KF_y = KF_y + KF_K * (KF_H * [obj.eye_records(i,obj.Y_COORD); KF_y_velocity_degree(i)] - KF_H * KF_y);
                    KF_P = KF_P - KF_K * KF_H * KF_P;
                end
            end
            
            for i=1:length(obj.eye_records)
                KF_window_start = i;
                KF_window_ended = i+obj.window_size;
                if(KF_window_ended > length(obj.eye_records) )
                    break;
                end
                KF_chi=0;
                for j=KF_window_start:KF_window_ended
                    if( obj.eye_records( j,obj.VALIDITY) == obj.DATA_VALID )
                        KF_chi = KF_chi + ( KF_measured(j) - KF_predicted(j) ).^2 / obj.deviation;
                    end
                end
                if( obj.eye_records( KF_window_ended,obj.VALIDITY ) == obj.DATA_VALID)
                    if (abs(KF_chi) < obj.chi_threshold)
% Classify points below Chi-square threshold as fixations
                        obj.eye_records(KF_window_ended,obj.MOV_TYPE) = obj.FIXATION_TYPE;
                    else
% Classify points above Chi-square threshold as saccades
                        obj.eye_records(KF_window_ended,obj.MOV_TYPE) = obj.SACCADE_TYPE;
% Due to our valocity calculations (left difference) we should include
% previous point to the saccade too
                        if( KF_window_ended > 1)
                            if( obj.eye_records( KF_window_ended - 1,obj.VALIDITY ) == obj.DATA_VALID )
                                obj.eye_records( KF_window_ended - 1,obj.MOV_TYPE ) = obj.SACCADE_TYPE;
                            end
                        end
                    end
                end
            end
         
% Special case
            obj.eye_records( ( isnan(obj.eye_records(:,obj.MOV_TYPE)) ),obj.MOV_TYPE ) = obj.NOISE_TYPE;
            
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Complete data classification with I-KF classifier in :',datestr(now),'\n'));
            end
        end

% Access interface to publicly accessible properties
        function set.chi_threshold(obj, value)
            obj.chi_threshold = value;
        end
        
        function set.window_size(obj, value)
            obj.window_size = value;
        end
        
        function set.deviation(obj, value)
            obj.deviation = value;
        end
    end
end
