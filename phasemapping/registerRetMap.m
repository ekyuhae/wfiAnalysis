function cord = registerRetMap(sPath, barmapping, expID_ret)

addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
idx = strfind(sPath, filesep);
expID= sPath(idx(end-1)+1:idx(end)-1);

% expID= sPath(idx(end-2)+1:idx(end-1)-1);
Animal = regexp(sPath, '(?<=\\WFI\\)M\d{6}_\d', 'match');
Animal = Animal{1};

if barmapping 
    rawPath = ['Y:\haider\Data\WFI\Mapping\Animals\' Animal filesep 'barmapping' filesep expID];
else 
    rawPath = ['Y:\haider\Data\WFI\Mapping\Animals\' Animal filesep 'gratings' filesep expID];
end 
snap = imread([rawPath filesep 'Snapshot_1.jpg']);

cd(['Y:\haider\Data\analyzedData\EKK\WFI\' Animal filesep 'retinotopy' filesep expID_ret])
load('plotPhaseMap.mat');
load('plotAmpMap.mat');
load('cPhaseMaps.mat');
snap_ret = load('Snapshot_1.mat');
snap_ret = snap_ret.snap; 

% snap_ret=imresize(snap_ret, [size(plotPhaseMap,1) size(plotPhaseMap,2)]);
% snap = imresize(snap, [size(plotPhaseMap,1) size(plotPhaseMap,2)]);

figure
imagesc(snap_ret);axis image; colormap gray; freezeColors;
hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([-0.5 0.5]);
set(vfsIm,'AlphaData',plotAmpMap);axis image;
set(gca, 'XTick', [], 'YTick', []);

ret=getframe;
% imwrite(ret.cdata, [sPath filesep 'vesselmap.png'])
axis off

ret = frame2im(ret);
[combined, ~,~, cord] = imalign_AR(imadjust(snap), ret);
save([sPath filesep 'coordinates.mat'], 'cord', 'snap')

VFS_Outline_Contour(plotPhaseMap, cPhaseMaps{2,1}, Animal, sPath);

% figure;
% imshow(combined)
