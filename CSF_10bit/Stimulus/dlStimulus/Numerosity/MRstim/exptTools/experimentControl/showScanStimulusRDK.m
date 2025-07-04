function [response, timing, quitProg] = showScanStimulusRDK(display,stimulus, t0)
% [response, timing] = showStimulus(display,stimulus, time0)
%
% t0 is the time the scan started and the stimulus timing should be
% relative to t0. If t0 does not exist it is created at the start of
% this program.
%  
% HISTORY:
% 2005.02.23 RFD: ported from showStimulus.
% 2005.06.15 SOD: modified for OSX. Use internal clock for timing rather
% than framesyncing because getting framerate does not always work. Using
% the internal clock will also allow some "catching up" if stimulus is
% delayed for whatever reason. Loading mex functions is slow, so this 
% should be done before callling this program.

% input checks
if nargin < 2,
	help(mfilename);
    return;
end;
if nargin < 3 || isempty(t0),
    t0 = GetSecs; % "time 0" to keep timing going
end;


% quit key
try 
    quitProgKey = display.quitProgKey;
catch
    quitProgKey = KbName('q');
end;

% some variables
nFrames = size(stimulus.images, 3);
HideCursor;
nGamma = size(stimulus.cmap,3);
response.keyCode = zeros(length(stimulus.seq),2); % get 1 buttons max
response.secs = zeros(size(stimulus.seq));        % timing
quitProg = 0;

rect=Screen('Rect', display.windowPtr);
if ~isfield(display, 'Rect');
    tmp1=round((display.numPixels(1)-display.numPixels(2))/2);
    display.Rect=[tmp1; 0; tmp1+display.numPixels(2); display.numPixels(2)];    
end

stimulus.makeMovie=1;
if isfield(stimulus, 'makeMovie') && stimulus.makeMovie==1
    aviobj = avifile('dotsMovie.avi', 'FPS', 20);
end

% go
disp(sprintf('[%s]:Running. Hit %s to quit.',mfilename,KbName(quitProgKey)));
for frame = 1:nFrames	
    
    %--- update display
    % If the sequence number is positive, draw the stimulus and the
    % fixation.  If the sequence number is negative, draw only the
    % fixation.

   % if stimulus.seq(frame)>0
        % put in an image 
		%imgNum = mod(stimulus.seq(frame)-1,nImages)+1;
        dotPattern=stimulus.images(:,:,frame);
        dotPattern=dotPattern(dotPattern<55555);
        dotPattern=reshape(dotPattern, 2, length(dotPattern)/2); 
        
        cols=stimulus.images(1,:,frame)<55555;
        dotColors=stimulus.seq.colors(cols==1);
        dotColors=repmat(dotColors,3,1);
        
        
        Screen('FillRect', display.windowPtr, [128 128 128], rect );
       
        
        if size(dotPattern,2)>1
            Screen('DrawDots',display.windowPtr, double(dotPattern), double(stimulus.seq.size), double(dotColors), [display.Rect(1) display.Rect(2)],1);
        end
        
        %Screen('DrawTexture', display.windowPtr, stimulus.textures(imgNum), stimulus.srcRect, stimulus.destRect);
        
        drawFixation(display,stimulus.fixSeq(frame));
%     elseif stimulus.seq(frame)<0
%         % put in a color table
% 		gammaNum = mod(-stimulus.seq(frame)-1,nGamma)+1;
%         % The second argument is the color index.  This apparently changed
%         % in recent times (07.14.2008). So, for now we set it to 1.  It may
%         % be that this hsould be 
%         drawFixation(display,stimulus.fixSeq(frame));
% 		Screen('LoadNormalizedGammaTable', display.windowPtr, stimulus.cmap(:,:,gammaNum));
%     end;
    
    %--- timing
    waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
    
    %--- get inputs (subject or experimentor)
    while(waitTime<0),
        % Scan the UMC device for subject response
        [ssKeyCode,ssSecs] = deviceUMC('response_and_trigger',display.devices.UMCport);
        if any(ssKeyCode)
            response.keyCode(frame,:) = ssKeyCode; 
            response.secs(frame)    = ssSecs - t0;
        end;
        % scan the keyboard for experimentor input
        [exKeyIsDown,~,exKeyCode] = KbCheck(display.devices.keyInputInternal);
        if(exKeyIsDown)
            if(exKeyCode(quitProgKey)),
                quitProg = 1;
                break; % out of while loop
            end;
        end;

        % if there is time release cpu
        if(waitTime<-0.02),
            WaitSecs(0.01);
        end;
        
        % timing
        waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
    end;
    
    %--- stop?
    if quitProg,
        disp(sprintf('[%s]:Quit signal recieved.',mfilename));
        break;
    end;

    %--- update screen
    Screen('Flip',display.windowPtr);
    
    if isfield(stimulus, 'makeMovie') && stimulus.makeMovie==1 && frame<=600;
        if isfield(display, 'Rect')
            imageArray=Screen('GetImage', display.windowPtr, display.Rect);
        else
            imageArray=Screen('GetImage', display.windowPtr);
        end
        aviobj = addframe(aviobj,imageArray);
    end
    if isfield(stimulus, 'makeMovie') && stimulus.makeMovie==1 && frame==601;
        aviobj = close(aviobj);
    end
end;

% that's it
ShowCursor;
timing = GetSecs-t0;
disp(sprintf('[%s]:Stimulus run time: %f seconds [should be: %f].',mfilename,timing,max(stimulus.seqtiming)));

return;
