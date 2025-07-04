
% ret - program to start retinotopic mapping experiments (under OSX)
%       
%
% 06/2005 SOD Ported to OSX. If the mouse is invisible,
%             moving it to the Dock usually makes it reappear.
% 10/2005 SOD Several changes, including adding gui.

% clean up - it's good to clean up but mex files are extremely slow to be
% loaded for the first time in MacOSX, so I chose not to do this to speed
% things up.
%close all;close hidden;
%clear mex;clear all;
%pack;

%set to 1 to use with ViewPoint eyetracker (Utretcht UMC fMRI eyetracker)
eyeTracker=0;

% ask parameters with user interface

params = retMenu;
drawnow;
params.NatVPhs = 2;
params = setRetinotopyParams(params.experiment, params);

% set response
% deviceyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
% old way keep for compatibility
params.devices = getDevices(false);


if isempty(params.devices.keyInputExternal)
    params.devices.keyInputExternal = params.devices.keyInputInternal;
end
% new way
params.devices.UMCport = -1; % default = 2; set to -1 for testing without device
if params.devices.UMCport ==2
    diary on
end

fprintf('[%s]:Using port %d.\n',mfilename,params.devices.UMCport);
deviceUMC('open',params.devices.UMCport);

if eyeTracker==1
    params.eyeTracker=1;
else
    params.eyeTracker = 0;
end

% go




oldLevel = Screen('Preference', 'Verbosity', 1);

doRetinotopyScan(params);
Screen('Preference', 'Verbosity', oldLevel);

% clean up
deviceUMC('close',params.devices.UMCport);
