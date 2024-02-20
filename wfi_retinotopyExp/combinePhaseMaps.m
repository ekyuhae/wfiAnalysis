for i = 1:2 
cPhaseMaps11_scaled{i,1} = imresize(cPhaseMaps11{i,1}, size(cPhaseMaps{i,1}));
combinedPM{i,1} = cat(3, cPhaseMaps{i,1}, cPhaseMaps11_scaled{i,1});
avgPM{i,1} = mean (combinedPM{i,1}, 3);
end 

cMagMaps11_scaled{1,1} = imresize(cMagMaps11{1,1}, size(cMagMaps{1,1}));
combinedMM{1,1} = cat(3, cMagMaps{1,1}, cMagMaps11_scaled{1,1});
avgMM{1,1} = mean (combinedMM{1,1}, 3);

[dhdx, dhdy] = gradient(avgPM{1,1}); % get gradient of horizontal phasemap
[dvdx, dvdy] = gradient(avgPM{2,1});

graddir_hor = atan2(dhdy,dhdx); % get gradient direction angle 
graddir_vert = atan2(dvdy,dvdx);
vdiff = exp(-1i*graddir_hor) .* exp(1i*graddir_vert); % changed to vertical -horizontal to get blue for V1 

VFS{1} = sin(angle(vdiff)); %Visual field sign map
VFS{1} = spatialFilterGaussian(VFS{1},2);
cVFS{1,1} = median(cat(3,VFS{:}),3);

h = figure;
subplot(2,2,1);
imagesc(avgPM{1,1});axis image; colormap hsv; colorbar; %freezeColors;
title('Elevation - nTrials = 33'); caxis([-40 40])
subplot(2,2,2);
imagesc(avgPM{2,1});axis image; colormap hsv; colorbar; %freezeColors;
title('Azimuth - nTrials = 33');  caxis([-15 120])
subplot(2,2,3);
imagesc(avgMM{1,1});axis image; colorbar; colormap jet;
title('Mean Magnitude');
subplot(2,2,4);
imagesc(spatialFilterGaussian(cVFS{1,1},2)); axis image;colorbar
caxis([-0.5 0.5]) % apply colorbar threshold
title(['VisualFieldSign - binSize = ' num2str(1) '; smth = ' num2str(phaseMapSmth)]);
savefig(h,[img.savePath filesep Animal '_phaseMap_allPlots_33_trials.fig']);
h.PaperUnits = 'inches';
set(h, 'PaperPosition', [0 0 15 15]);
saveas(h,[img.savePath filesep Animal '_phaseMap_allPlots_33_trials.jpg'])
clear h

VFS_Outline(spatialFilterGaussian(cVFS{1,1},2),avgPM{2,1},avgPM{1,1}, Animal, img.savePath);

    %nested functions
function img = spatialFilterGaussian(img, sigma)
if sigma > 0 && (numel(img) ~=  sum(sum(isnan(img))))
    hh = fspecial('gaussian',size(img),sigma);
    hh = hh/sum(hh(:));
    img = ifft2(fft2(img).*abs(fft2(hh)));
end
end