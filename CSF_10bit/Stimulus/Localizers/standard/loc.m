% loc - UPDATE MARCUS
%       
% 06/2005 SOD Ported to OSX. If the mouse is invisible,
%             moving it to the Dock usually makes it reappear.
% 10/2005 SOD Several changes, including adding gui.
% 04/2006 SOD Converted from ret to loc

% clean up - it's good to clean up but mex files are extremely slow to be
% loaded for the first time in MacOSX, so I chose not to do this to speed
% things up. `
% close all;close hidden;
% clear mex;clear all;
close all
clear all
sca;
Screen('Preference', 'SkipSyncTests', 1); 
CSF10_bit_path = '/data1/projects/dumoulinlab/Lab_members/Marcus/programs/Experiments/vdncsf_protocol/CSF_10bit';
addpath(genpath(CSF10_bit_path))

% add stimulus files

commandwindow
% pack; 

% ask parameters with user interface
params = locMenu;
drawnow;

% now set rest of the params
params = setLocParams(params.experiment, params);

% set response device
params.devices = getDevices(false);

if isempty(params.devices.keyInputExternal)
    params.devices.keyInputExternal = params.devices.keyInputInternal;
end
fprintf('[%s]:Getting subjects responses from device #%d\n',mfilename,params.devices.keyInputExternal);
fprintf('[%s]:Getting experimentor''s responses from device #%d\n',mfilename,params.devices.keyInputInternal);

% go
doLocScan(params);
sca;

