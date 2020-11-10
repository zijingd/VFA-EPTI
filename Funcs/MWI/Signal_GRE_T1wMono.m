%% S = Signal_GRE_T1wMono(Mo, FA, T1, TR,varargin)
%
% Input:
% ------
%	Mo: Proton density, can be arbitary unit
%	FA: Flip angle
% 	T1: (same unit as TR)
%	TR: (same unit as T1)
% Output:
% -------
%	S: T1-weighted signal
%
% Description: Calculation of signal intensity with known flip angle, T1 and TR
%
%   Author: Kwok-shing Chan @ University of Aberdeen
%   Date created: Jan 1, 2016
%   Ref: Rapid combined T1 and T2 mapping using gradient recalled
%   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%   Date last edited: 20 February 2018
%
function S = Signal_GRE_T1wMono(rho,FA,T1,TR,varargin)
% parse input arguments
[b1,displayImage,displaySurf,displayPlot,verbose,isDictionary]=parse_varargin_Signal_GRE_T1wMono(varargin);

if length(T1)>1 && length(FA)>1 && verbose
    disp('More than one T1. Output will be in 2D.');
end
if length(rho)<length(T1)
    disp('No. of rhoes is less than no of T1s, assuming all T1s have the same rho');
    rho = rho*ones(length(T1));
end

% if B1 is provided, apply it here
FA = FA.*b1;

%% Core algorithm
if ~isDictionary
    S = zeros(length(T1),length(FA));
    for k=1:length(T1)
        E1 = exp(-TR/T1(k));
        S(k,:) = (rho(k).*(1-E1).*sind(FA))./(1-E1.*cosd(FA));
    end
else
    E1 = exp(-TR./T1);
    S = (rho.*(1-E1).*sind(FA))./(1-E1.*cosd(FA));
end

%% for simulation display
if displayImage
    figure;imagesc(S);xlabel('Flip angle (degree)');ylabel('T1');colormap jet;
    set(gca,'Ydir','normal');
    yTickLabelName = get(gca,'yticklabel');
    for kk = 1:length(yTickLabelName)
        yTickLabelName{kk,1} = num2str(T1(str2double(yTickLabelName{kk,1})));
    end
    set(gca,'yTickLabel',yTickLabelName);
    xTickLabelName = get(gca,'xticklabel');
    for kk = 1:length(xTickLabelName)
        xTickLabelName{kk,1} = num2str(FA(str2double(xTickLabelName{kk,1})));
    end
    set(gca,'xTickLabel',xTickLabelName);
    title(sprintf('SPGR T1 weighting, TR=%f',TR));
end
if displaySurf
    [rx,ry] = meshgrid(FA,T1);
    figure;surf(rx,ry,S);xlabel('Flip angle (degree)');ylabel('T1');colormap jet;
    title(sprintf('SPGR T1 weighting, TR=%f',TR));
end
if displayPlot
    figure;
    for k = 1:length(T1)
        plot(FA,squeeze(S(k,:))/sum(rho(:)));hold on;
    end
    xlabel('Flip angle (degree)');ylabel('% intensity change to total signal');colormap jet;
    title(sprintf('SPGR T1 weighting, TR=%f',TR));
end
end

function [b1,displayImage,displaySurf,displayPlot,verbose,isDictionary]=parse_varargin_Signal_GRE_T1wMono(arg)
% predefine parameters
displayImage = false;    % in second
displaySurf = false;
displayPlot = false;
verbose = false;
isDictionary = false;
b1 = 1;

% look for flags and 'Name/Value' pairs
for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'image')
        displayImage = true;
    end
    if strcmpi(arg{kvar},'surf')
        displaySurf = true;
    end
    if strcmpi(arg{kvar},'plot')
        displayPlot = true;
    end
    if strcmpi(arg{kvar},'-v')
        verbose = true;
    end
    if strcmpi(arg{kvar},'dict')
        isDictionary = true;
    end
    if strcmpi(arg{kvar},'b1')
        b1 = arg{kvar+1};
    end
end
end