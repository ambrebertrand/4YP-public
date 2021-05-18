function [output,values_for_csv] = histograms(roi_vals_vis1, roi_vals_vis2, modality_name, modality_units)

%% script to extract voxel ROI data from Mirada output
% Daniel Bulte, IBME, University of Oxford, December 2020

% edited by Ambre Bertrand, February 2021
% use this version as part of the registration pipeline


%% set vector to plot fits
rango = transpose(linspace(600,1600,200));

%% ROI for visit 1
roivals1 = double(roi_vals_vis1(:));
roivals1 = roivals1(roivals1~=0);
roivals1(le(roivals1,0))=nan;
numvox1 = length(roivals1(:));

% fit components GMMs to data and select best one
k_max = 3;   % 3 tends to work ok

AIC_1 = zeros(1,k_max);
GMModels_1 = cell(1,k_max);
options = statset('MaxIter',500);

for k = 1:k_max
    GMModels_1{k} = fitgmdist(roivals1,k,'Options',options,'Regularization',0.2);
    AIC_1(k)= GMModels_1{k}.AIC;
end

figure
plot(AIC_1)
title(sprintf('AIC scores against number of GMM components \n for %s tumour data Visit 1', modality_name));
xlabel('k')
ylabel('AIC score')
xticks(1:k_max)


[minAIC_1,numComponents_1] = min(AIC_1);

GMModel_1 = GMModels_1{numComponents_1};

ExpCovariance_1(1:numComponents_1) = GMModel_1.Sigma(:,:,1:numComponents_1);
ExpSigma_1 = sqrt(ExpCovariance_1);


% find mean and stddev of all values
mean1 = mean(roivals1,'omitnan');
std1 = std(roivals1,'omitnan');
skew1 = skewness(roivals1);
kurt1 = kurtosis(roivals1);

gm1 = gmdistribution(GMModel_1.mu,GMModel_1.Sigma,GMModel_1.ComponentProportion);
gmPDF1 = pdf(gm1,rango);

figure;
h1 = histogram(roivals1,50);
hold on
title(sprintf('ROI values for %s Visit 1, GMM with %d components', modality_name, numComponents_1))
plot(rango,h1.BinWidth*numvox1*gmPDF1,'LineWidth',2.0);



%% repeat above for second visit

roivals2 = double(roi_vals_vis2(:));
roivals2 = roivals2(roivals2~=0);
roivals2(le(roivals2,0))=nan;
numvox2 = length(roivals2(:));

% fit 1, 2 & 3 component GMM's to data and select best one
k_max = 3; 

AIC_2 = zeros(1,k_max);
GMModels_2 = cell(1,k_max);
options = statset('MaxIter',500);

for k = 1:k_max
    GMModels_2{k} = fitgmdist(roivals2,k,'Options',options,'Regularization',0.2);
    AIC_2(k)= GMModels_2{k}.AIC;
end

figure
plot(AIC_2)
title(sprintf('AIC scores against number of GMM components \n for %s tumour data Visit 2', modality_name));
xlabel('k')
ylabel('AIC score')
xticks(1:k_max)


[minAIC_2,numComponents_2] = min(AIC_2);

GMModel_2 = GMModels_2{numComponents_2};

ExpCovariance_2(1:numComponents_2) = GMModel_2.Sigma(:,:,1:numComponents_2);
ExpSigma_2 = sqrt(ExpCovariance_2);

mean2 = mean(roivals2,'omitnan');
std2 = std(roivals2,'omitnan');
skew2 = skewness(roivals2);
kurt2 = kurtosis(roivals2);

gm2 = gmdistribution(GMModel_2.mu,GMModel_2.Sigma,GMModel_2.ComponentProportion);
gmPDF2 = pdf(gm2,rango);

figure;
h2 = histogram(roivals2,50);
hold on
title(sprintf('ROI values for %s Visit 2, GMM with %d components', modality_name, numComponents_2))
plot(rango,h2.BinWidth*numvox2*gmPDF2,'LineWidth',2.0);


%% plot both histograms on one figure for comparison
figure
h1 = histogram(roivals1,50);
hold on
h1=plot(rango,h1.BinWidth*numvox1*gmPDF1,'LineWidth',2.0);
hold on
h2 = histogram(roivals2,50);
hold on
title(sprintf('ROI values for %s, Visits 1 and 2', modality_name))
h2=plot(rango,h2.BinWidth*numvox2*gmPDF2,'LineWidth',2.0);
legend([h1,h2],{'Visit 1','Visit 2'});
xlabel([modality_units]);
ylabel('Pixel frequency');


%% Determine degree of response by reduction of low ADC peaks

outprop1 = zeros(1,3);
outmu1 = zeros(1,3);
outsig1 = zeros(1,3);

lowvol_1 = 0;
for comp = 1:numComponents_1
    if GMModel_1.mu(comp)<1200 % value was chosen arbitrarily and needs to be determined properly
        lowvol_1 = lowvol_1 + numvox1*GMModel_1.ComponentProportion(comp);
    end
    outprop1(comp) = GMModel_1.ComponentProportion(comp);
    outmu1(comp) = transpose(GMModel_1.mu(comp));
    outsig1(comp) = ExpSigma_1(comp);
end

outprop2 = zeros(1,3);
outmu2 = zeros(1,3);
outsig2 = zeros(1,3);

lowvol_2 = 0;
for comp = 1:numComponents_2
    if GMModel_2.mu(comp)<1200 % value was chosen arbitrarily and needs to be determined properly
        lowvol_2 = lowvol_2 + numvox2*GMModel_2.ComponentProportion(comp);
    end
    outprop2(comp) = GMModel_2.ComponentProportion(comp);
    outmu2(comp) = transpose(GMModel_2.mu(comp));
    outsig2(comp) = ExpSigma_2(comp);
end

deltalowpeak = 100*(lowvol_2-lowvol_1)./lowvol_1; 

deltamean = 100*(mean2-mean1)./mean1;

deltavol = 100*(numvox2-numvox1)./numvox1;

output = [numvox1 mean1 std1 outprop1 outmu1 outsig1 numvox2 mean2 std2 outprop2 outmu2 outsig2 deltavol deltamean deltalowpeak];

% 'Visit','Tumour volume (cm^3)','Biomarker','Mean','SD','Number of GMM components','Component 1 mean','Component 2 mean','Component 3 mean','Component 1 prop','Component 2 prop','Component 3 prop'
values_for_csv = [1,0,0,mean1,std1,numComponents_1,outmu1(1),outmu1(2),outmu1(3),outprop1(1),outprop1(2),outprop1(3);
                  2,0,0,mean2,std2,numComponents_2,outmu2(1),outmu2(2),outmu2(3),outprop2(1),outprop2(2),outprop2(3)];
% keep 2nd and 3rd index empty for tumour volume and modality name
              
end 


% FOR EACH MODALITY at visit 1 and visit 2
% number of gmm components
% prop of components
% mean of each components
% general mean
% std 
% tumour volume


