% add stimulus files
CSF10_bit_path = '/data1/projects/dumoulinlab/Lab_members/Marcus/programs/Experiments/vdncsf_protocol/CSF_10bit';
% addpath(genpath(CSF10_bit_path))
CSF10_bit_ptb_path = '/packages/matlab/toolbox/psychtoolbox/3.0.17';
% add psychtoolbox 
addpath(genpath(CSF10_bit_ptb_path));
cd(CSF10_bit_ptb_path);
SetupPsychtoolbox;
cd(CSF10_bit_path);

%% Setup 14 bit screen, 7T
Screen('Preference', 'SkipSyncTests', 1); % ADDED BY MARCUS - TESTING
% press 'b' when asked to input a letter.
BitsPlusImagingPipelineTest(2)

BitsPlusIdentityClutTest(1, 0)
% Wait for a bit. 
% Press space until the stimulus on screen matches the description on
% screen.
% Then press 's', ' Esc', 'y', 'Enter' to save and confirm the screen settings 
% BitsPlusIdentityClutTest(1, 0)

%% Setup 14 bit screen, eyetracker room
% press 'b' when asked to input a letter.
% BitsPlusImagingPipelineTest(1)
% 
% BitsPlusIdentityClutTest(1, 0)
% Wait for a bit. 
% Press space until the stimulus on screen matches the description on
% screen.
% Then press 's', ' Esc', 'y', 'Enter' to save and confirm the screen settings 

