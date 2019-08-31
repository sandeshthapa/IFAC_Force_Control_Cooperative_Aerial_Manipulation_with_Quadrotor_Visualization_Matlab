clc;
clear all;
close all;

%   obj = VideoReader('myVideo2.avi');
%     vid = readFrame(obj);
%     nFrame = size(vid,122);  
%     for k = 1 : nFrame
%         newname = strcat(num2str(k),'.jpg');
%         imwrite(vid(:,:,:,k), newname);  
%         pause(0.1);
%     end  

    
 vid=VideoReader('myVideo2.avi');
 numFrames = vid.NumberOfFrames;
 n=numFrames;
 for i = 1:1:n
 frames = read(vid,i);
 imwrite(frames,['Frame\Image' int2str(i), '.png']);
 im(i)=image(frames);
 end

% vidObj = VideoReader('myVideo2.avi');
% vidObj.CurrentTime = 0.1;  
% frames = vidObj.CurrentTime; 
% for f = 1: frames
%     thisframe = readFrame(vidObj,f);
%     figure(1); 
%     imagesc(thisframe);
%     thisfile = sprintf('frame_%04d.jpg',f);
%     imwrite(thisframe,thisfile);
% end 

% 
% vidHeight = vidObj.Height; 
% vidWidth = vidObj.Width; 
% 
% s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
%     'colormap',[]);
% 
% s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
%     'colormap',[]);
% 
% whos s
% image(s(1).cdata)

% set(gcf,'position',[150 150 vidObj.Width vidObj.Height]);
% set(gca,'units','pixels');
% set(gca,'position',[0 0 vidObj.Width vidObj.Height]);
% movie(s,1,vidObj.FrameRate);
% 


% vidObj.CurrentTime = 0.1; 
% currAxes = axes;
% while hasFrame(vidObj)
%     vidFrame = readFrame(vidObj);
% %     image(vidFrame, 'Parent', currAxes);
% %     currAxes.Visible = 'off';
% %     pause(1/v.FrameRate);
% end
% 
% whos vidFrame