function doLocScan(params)
% doLocScan(params)
%
% Runs any of several retinotopy scans
%
% 99.08.12 RFD wrote it, consolidating several variants of retinotopy scan code.
% 05.06.09 SOD modified for OSX, lots of changes.
% 04/2006  SOD converted from doRetinotopyScan

% defaults
if ~exist('params', 'var')
	error('No parameters specified!');
end

% quit key
try
    quitProgKey = params.display.quitProgKey;
catch
    quitProgKey = KbName('q');
end

% make/load stimulus
switch params.experiment
    case 'moving bars vs fixation'
        [stimulus, onebackSequence] = makeBarStimulus(params);
    case {'contrast sensitivity function 4sf', 'contrast sensitivity function 6sf'}
        [stimulus, onebackSequence, params] = makeCrfStimulus(params);
    otherwise
        [stimulus, onebackSequence] = makeLocStimulus(params);
end

% add scotoma
params.scotoma = []; % [1.46 -4.91 0.565 0.361 131 5];
if isfield(params,'scotoma')
    if ~isempty(params.scotoma)
        stimulus = insertScotoma(stimulus,params,params.scotoma);
    end
end

% loading mex functions for the first time can be
% extremely slow (seconds!), so we want to make sure that 
% the ones we are using are loaded.
KbCheck;GetSecs;WaitSecs(0.001);

try
    % check for OpenGL
    AssertOpenGL;
    
    Screen('Preference','SkipSyncTests', 1);
    
    % Open the screen
    params.display                = openScreen(params.display);
    params.display.devices        = params.devices;

    % to allow blending
    % Screen('BlendFunction', params.display.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Store the images in textures
    stimulus = createTextures(params.display,stimulus);

    % getting to the source rather than doScan->doTrial->showStimulus 
    for n = 1:params.repetitions
        % set priority
%        Priority(params.runPriority);
        
        % wait for go signal
        % waituntillspacepress;
        % pressKey2Begin(params.display);      

        % countdown + get start time (time0)
        % [time0] = countDown(params.display,params.startScan+2,params.startScan);

        % Trigger ret code
%        dispStringInCenter(params.display,'Waiting for scanner to start',0.6,[],'white');
        
        dispStringInCenter(params.display,'Waiting for scanner to start',0.6,40,[0 0 0]);
        
        KbQueueRelease() 
        KbTriggerWait(KbName('t')); %,-3)%-3 makes sure the Forb box is caught
        % [junk time0]=deviceUMC('wait for trigger',params.devices.UMCport);  
   
        % go
        %[response, timing, quitProg] = showScanStimulus(params.display,stimulus,time0);
        [response, timing, quitProg] = showScanStimulus(params.display,stimulus);
        
        % reset priority
        % Priority(0);
        
        % get performance
        [pc,rc,nn] = getFixationPerformance(params.fix,stimulus,response); 
        % disp(sprintf('[%s]:Fixation dot task(%d): percent correct: %.1f %%, reaction time: %.1f secs',mfilename,nn,pc,rc));
        fprintf('[%s]:Fixation dot task(%d): percent correct: %.1f %%, reaction time: %.1f secs',mfilename,nn,pc,rc);
        
%         % get 1-back performance
%         if exist('onebackSequence') && ~isempty(onebackSequence)
%             [pc,rc,nn] = getDetectionPerformance(params.fix,stimulus,response,onebackSequence);
%             % disp(sprintf('[%s]:One-back task(%d): percent correct: %.1f %%, reaction time: %.1f secs',mfilename,nn,pc,rc));
%             fprintf('[%s]:One-back task(%d): percent correct: %.1f %%, reaction time: %.1f secs',mfilename,nn,pc,rc);
%         end
%         
        % save 
        if params.savestimparams == 1

            filename = [ '/data1/projects/dumoulinlab/Lab_members/Marcus/programs/Experiments/vdncsf_protocol/CSF_10bit/Output/' datestr(now,30) '.mat'];
            save(filename);                % save parameters
            % disp(sprintf('[%s]:Saving in %s.',mfilename,filename));
            fprintf('[%s]:Saving in %s.',mfilename,filename);
        end
        
        % keep going?
        if quitProg % don't keep going if quit signal is given
            break;
        end
    end
    
    % Close the one on-screen and many off-screen windows
    closeScreen(params.display);
catch
    % clean up if error occurred
    Screen('CloseAll');
%    setGamma(0);
    Priority(0);
    ShowCursor;
    rethrow(lasterror);
end
