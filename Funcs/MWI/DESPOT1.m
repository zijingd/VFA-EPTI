%% [t1, m0] = DESPOT1(S,FA,TR,varargin)
% 
% Input
% ----------
% S         : Set of T1-weighted signals from different flip angles, pixel-based
% FA        : Set of flip angles used corresponding to S
% TR        : TR of sequence (s)
% Flags     :
%   'range'     -   boundary of T1
%   'b1'        -   B1 map
%   'option'    -   fitting option
%   **solver choices**
%   'regression'	: using mldivide to solve b=x\y (default) (fastest)
%   'lsqnonneg'     : nonnegative linear least-squares (fast)
%   'lsqnonlin'     : nonlinear least square (slow)
%
% Output
% ----------
% t1        : Calcuted T1 (same as TR)
% m0        : Calculated proton density
%
% DESPOT1 to estimate T1 and proton density Mo
% Linear fitting of y = ax+b
% Rapid combined T1 and T2 mapping using gradient recalled acquisition in the steady state.
% Deoni et al. MRM 2003;49:515-526
%
% Kwok-shing Chan @ dccn
% k.chan@donders.ru.nl
% Date created: Jan 11, 2016
% Date last edited: 20 February, 2017
%
function [t1, m0] = DESPOT1(S,FA,TR,varargin)
% check error
if length(S) < 2
    error('Only one point to fit.');
end
if length(S) ~= length(FA)
    error('Size of signal is not equal to size of flip angle');
end

% parse argument input
[solver,b1,rangeT1,option] = parse_varargin_DESPOT1(varargin);

%% Core
% DESPOT1 formulation
S = double(S);
FA = double(FA)*b1;
TR = double(TR);

y = S(:)./sind(FA(:));
xCol = S(:)./tand(FA(:));
x = ones(length(S),2);
x(:,1) = xCol;
    
% solve DESPOT1 using different approach
switch solver
    case 'regression'
        b = x\y;
        t1 = -TR/log(b(1));
        m0 = b(2)/(1-exp(-TR/t1));
    case 'lsqnonneg'
        b = lsqnonneg(x,y);
        t1 = -TR/log(b(1));
        m0 = b(2)/(1-exp(-TR/t1));
    case 'lsqnonlin'
        % Obtain initial guesses
        b = x\y;
        t10 = -TR/log(b(1));
        m00 = b(2)/(1-exp(-TR/t10));
        T1_lb = rangeT1(1);
        T1_ub = rangeT1(2);
        c0 = [t10, m00];
        lb = [T1_lb, min(S)];
        ub = [T1_ub, 2*m00];

        res = lsqnonlin(@(x)fitError_DESPOT1(x,S,FA,TR),c0,lb,ub,option);
        t1 = res(1);
        m0 = res(2);
end
t1 = real(t1);
m0 = real(m0);
end

%% lsqnonlin: compute fitting residual 
function [fiter] = fitError_DESPOT1(x,S_meas,FA,TR)
% grab estimates
T1 = x(1);
Mo = x(2);

% simulate signal
S_fit = Signal_GRE_T1wMono(Mo,FA,T1,TR);
% compute fiter, using magnitude fitting
fiter = computeFiter(S_meas,S_fit,length(S_fit));

end

%% parse argument input
function [solver,b1,rangeT1,option] = parse_varargin_DESPOT1(arg)
% predefine parameters
rangeT1 = [0,5];    % in second
b1 = 1;
option = [];
solver = 'regression';

% look for flags and 'Name/Value' pairs
for kvar = 1:length(arg)
%     if strcmpi(arg{kvar},'regression')
%         solver = 'regression';
%     end
    if strcmpi(arg{kvar},'lsqnonlin')
        solver = 'lsqnonlin';
    end
    if strcmpi(arg{kvar},'lsqnonneg')
        solver = 'lsqnonneg';
    end
    if strcmpi(arg{kvar},'b1')
        b1 = arg{kvar+1};
    end
    if strcmpi(arg{kvar},'range')
        rangeT1 = arg{kvar+1};
    end
    if strcmpi(arg{kvar},'option')
        option = arg{kvar+1};
    end
end
% lsqnonlin fitting option
if isempty(option) && strcmpi(solver,'lsqnonlin')
    option = optimoptions(@lsqnonlin,'Display','off','Jacobian','off','DerivativeCheck','off','MaxIter',500);
end
end
