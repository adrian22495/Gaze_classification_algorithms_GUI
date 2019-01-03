% This script is just for testing purposes only - I don't know how to use
% test-unit or it's analog in MatLab. At least now, i hope. So I just use
% this mock script to check all my classes.
clear 
%%%%%%%%%%%%%%%%%%% IKF END 
IKF = classificator_IKF_class;
IKF.input_data_name                = 'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42.txt';
IKF.x_field =                      8;
IKF.y_field =                      9;
IKF.header_count =                 11;
IKF.read_data();
if(IKF.error_code ~= 0)
    errordlg(IKF.error_message,'File reader error.');
else
    IKF.image_width_mm =               400;
    IKF.image_height_mm =              300;
    IKF.image_width_etu =              1024;
    IKF.image_height_etu=              768;
    IKF.distance_from_screen =         800;
    IKF.chin_rest_height =             300;
    IKF.chin_rest_position =           200;
    IKF.convert_from_ETU_to_degrees();
    IKF.convert_from_degree_to_ETU();
    IKF.delta_t_sec =                  .001;
    IKF.chi_threshold =                 3.75;
    IKF.window_size =                   5;
    IKF.deviation =                     1000;
    IKF.classify();

    IKF.merge_records();

    IKF.minimal_saccade_amplitude =     1;
    IKF.maximal_saccade_amplitude =     24;
    IKF.minimal_saccade_length =        4;
    IKF.unfiltered_saccade_records = IKF.saccade_records;
    IKF.saccade_filtering();
    IKF.saccade_records = IKF.filtered_saccade_records;

    IKF.basename_output_filename =      'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42_ikf';
    IKF.basename_output_extension =     '.txt';
    IKF.setup_output_names();
    IKF.write_datafiles();

    fprintf('IKF fixations = %f saccade = %f pursuit = %f noise = %f\n',IKF.fixation_percentage,IKF.saccade_percentage,IKF.pursuit_percentage,IKF.noise_percentage);
end
clear IKF;
%%%%%%%%%%%%%%%%%%% IKF END 

%%%%%%%%%%%%%%%%%%% IVT BEGIN
IVT = classificator_IVT_class;
IVT.input_data_name                = 'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42.txt';
IVT.x_field =                       8;
IVT.y_field =                       9;
IVT.header_count =                  11;
IVT.read_data();
if(IVT.error_code ~= 0)
    errordlg(IVT.error_message,'File reader error.');
else
    IVT.image_width_mm =                400;
    IVT.image_height_mm =               300;
    IVT.image_width_etu =               1024;
    IVT.image_height_etu=               768;
    IVT.distance_from_screen =          800;
    IVT.chin_rest_height =              300;
    IVT.chin_rest_position =            200;
    IVT.convert_from_ETU_to_degrees();
    IVT.convert_from_degree_to_ETU();
    IVT.delta_t_sec =                   0.001;
    IVT.saccade_detection_threshold =   30;
    IVT.classify();

    IVT.merge_records();

    IVT.minimal_saccade_amplitude =     1;
    IVT.maximal_saccade_amplitude =     24;
    IVT.minimal_saccade_length =        4;
    IVT.unfiltered_saccade_records = IVT.saccade_records;
    IVT.saccade_filtering();

    IVT.saccade_records = IVT.filtered_saccade_records;

    IVT.basename_output_filename =      'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42_ivt';
    IVT.basename_output_extension =     '.txt';
    IVT.setup_output_names();
    IVT.write_datafiles();

    fprintf('IVT fixations = %f saccade = %f pursuit = %f noise = %f\n',IVT.fixation_percentage,IVT.saccade_percentage,IVT.pursuit_percentage,IVT.noise_percentage);
end
clear IVT;
%%%%%%%%%%%%%%%%%%% IVT END 

%%%%%%%%%%%%%%%%%%% IDT END 
IDT = classificator_IDT_class;
IDT.input_data_name                = 'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42.txt';
IDT.x_field =                      8;
IDT.y_field =                      9;
IDT.header_count =                 11;
IDT.read_data();
if(IDT.error_code ~= 0)
    errordlg(IDT.error_message,'File reader error.');
else
    IDT.image_width_mm =               400;
    IDT.image_height_mm =              300;
    IDT.image_width_etu =              1024;
    IDT.image_height_etu=              768;
    IDT.distance_from_screen =         800;
    IDT.chin_rest_height =             300;
    IDT.chin_rest_position =           200;
    IDT.convert_from_ETU_to_degrees();
    IDT.convert_from_degree_to_ETU();
    IDT.delta_t_sec =                  .001;
    IDT.dispersion_duration_sec_threshold = 0.1;
    IDT.dispersion_threshold =         0.6;
    IDT.classify();

    IDT.merge_records();

    IDT.minimal_saccade_amplitude =     1;
    IDT.maximal_saccade_amplitude =     24;
    IDT.minimal_saccade_length =        4;
    IDT.unfiltered_saccade_records = IDT.saccade_records;
    IDT.saccade_filtering();
    IDT.saccade_records = IDT.filtered_saccade_records;

    IDT.basename_output_filename =      'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42_idt';
    IDT.basename_output_extension =     '.txt';
    IDT.setup_output_names();
    IDT.write_datafiles();

    fprintf('IDT fixations = %f saccade = %f pursuit = %f noise = %f\n',IDT.fixation_percentage,IDT.saccade_percentage,IDT.pursuit_percentage,IDT.noise_percentage);
end
clear IDT;
%%%%%%%%%%%%%%%%%%% IDT END 

%%%%%%%%%%%%%%%%%%% IMST BEGIN
IMST = classificator_IMST_class;
IMST.input_data_name                = 'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42.txt';
IMST.x_field =                      8;
IMST.y_field =                      9;
IMST.header_count =                 11;
IMST.read_data();
if(IMST.error_code ~= 0)
    errordlg(IMST.error_message,'File reader error.');
else
    IMST.image_width_mm =               400;
    IMST.image_height_mm =              300;
    IMST.image_width_etu =              1024;
    IMST.image_height_etu=              768;
    IMST.distance_from_screen =         800;
    IMST.chin_rest_height =             300;
    IMST.chin_rest_position =           200;
    IMST.convert_from_ETU_to_degrees();
    IMST.convert_from_degree_to_ETU();
    IMST.delta_t_sec =                  .001;
    IMST.saccade_detection_threshold =  0.1;
    IMST.window_size =                   50;
    IMST.classify();

    IMST.merge_records();

    IMST.minimal_saccade_amplitude =     1;
    IMST.maximal_saccade_amplitude =     24;
    IMST.minimal_saccade_length =        4;
    IMST.unfiltered_saccade_records = IMST.saccade_records;
    IMST.saccade_filtering();
    IMST.saccade_records = IMST.filtered_saccade_records;

    IMST.basename_output_filename =      'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42_imst';
    IMST.basename_output_extension =     '.txt';
    IMST.setup_output_names();
    IMST.write_datafiles();
    
    fprintf('IMST fixations = %f saccade = %f pursuit = %f noise = %f\n',IMST.fixation_percentage,IMST.saccade_percentage,IMST.pursuit_percentage,IMST.noise_percentage);
end
clear IMST;
%%%%%%%%%%%%%%%%%%% IMST END 

%%%%%%%%%%%%%%%%%%% IHMM END 
IHMM = classificator_IHMM_class;
IHMM.input_data_name                = 'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42.txt';
IHMM.x_field =                      8;
IHMM.y_field =                      9;
IHMM.header_count =                 11;
IHMM.read_data();
if(IHMM.error_code ~= 0)
    errordlg(IHMM.error_message,'File reader error.');
else
    IHMM.image_width_mm =               400;
    IHMM.image_height_mm =              300;
    IHMM.image_width_etu =              1024;
    IHMM.image_height_etu=              768;
    IHMM.distance_from_screen =         800;
    IHMM.chin_rest_height =             300;
    IHMM.chin_rest_position =           200;
    IHMM.convert_from_ETU_to_degrees();
    IHMM.convert_from_degree_to_ETU();
    IHMM.delta_t_sec =                  .001;
    IHMM.saccade_detection_threshold =  30;
    IHMM.viterbi_sample_size =          100;
    IHMM.baum_welch_reiterations_count =4;
    IHMM.classify();

    IHMM.merge_records();

    IHMM.minimal_saccade_amplitude =     1;
    IHMM.maximal_saccade_amplitude =     24;
    IHMM.minimal_saccade_length =        4;
    IHMM.unfiltered_saccade_records = IHMM.saccade_records;
    IHMM.saccade_filtering();
    IHMM.saccade_records = IHMM.filtered_saccade_records;

    IHMM.basename_output_filename =      'C:\Users\EU\Desktop\Karpov\project\Input for scores\1D_HR_Saccade_participant-42_ihmm';
    IHMM.basename_output_extension =     '.txt';
    IHMM.setup_output_names();
    IHMM.write_datafiles();

    fprintf('IHMM fixations = %f saccade = %f pursuit = %f noise = %f\n',IHMM.fixation_percentage,IHMM.saccade_percentage,IHMM.pursuit_percentage,IHMM.noise_percentage);
end
clear IHMM;
%%%%%%%%%%%%%%%%%%% IHMM END 