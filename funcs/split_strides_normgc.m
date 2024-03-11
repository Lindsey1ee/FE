function strides = split_strides_normgc(variables, GE, npts, ploton)

% Splits a long time series of gait data into strides based on the gait
% events passed to the function. Also time normalizes each stride to have
% the same number of points so time can be plotted as a percentage of 
% the gait cycle.
%
% INPUTS:
%
% variables is a table or struct with the data to be split into strides
% GE are the gait events where each row is a stride with 5 gait events
%   RHS LTO LHS RTO RHSn
% npts is the number of points each normalized stride should have
% ploton is a flag about whether to plot the data. default is 1, to plot.

% OUTPUT
% 
% strides is a structure stride info/data
% if ploton = 1, then a figure for each variable will be generated

% checks number of input arguments. If only 3 variables passed, then set
% ploton to 1 to plot results
if nargin < 4, ploton=1; end 

% get names of variables passed to be split and time normalized
if isstruct(variables)
    varnames = fieldnames(variables);
elseif istable(variables)
    varnames = variables.Properties.VariableNames; 
end

strides.num = size(GE,1);
strides.lengths = GE.RHSn-GE.RHS+1;

%% loop through varnames and strides to split and time normalize
for v = 1:length(varnames)

    if ploton, figure, end

    for i = 1:strides.num

        % only run for "goodstrides." get_gaitevents puts 0's for gait
        % events for badstrides
        if GE.RHS(i) ~= 0
            % ~~ split data into strides ~~
            % get data for stride i and store into a cell array in a structure
            % field with the name of the variable. This was we can access the data later
            % need to use a cell array because each stride has different
            % lengths.
            strides.data.(varnames{v}){i} = variables.(varnames{v})(GE.RHS(i):GE.RHSn(i),:);

            % ~~ normalize lengths to npts ~~
            % all strides moving forward will have npts
            % use linspace to create a vector, x.
            % see linspace help to see what it does
            xi = linspace(1,strides.lengths(i),npts);

            % ~~ interpolate stride data at each value of xi ~~
            % see interp1 help to see what it does
            x = 1:strides.lengths(i);
            strides.normdata.(varnames{v})(:,i) = interp1(x,strides.data.(varnames{v}){i},xi);

            % ~~ plots individual strides ~~
            % see how each stride ends at different points
            if ploton
                subplot(1,2,1)
                hold on
                plot(strides.data.(varnames{v}){i})
                box off
                if i==1, xlabel("points"), ylabel(varnames{v}), end
            end
        end

    end

    % ~~ convert normadata to table if it is a structure ~~
    if isstruct(strides.normdata)
        strides.normdata = struct2table(strides.normdata);
    end

    % ~~ plot time normalized strides that all have the same points ~~
    if ploton
        strides.normgc = linspace(1,100,npts);
        subplot(1,2,2)
        plot(strides.normgc, strides.normdata.(varnames{v}))
        box off
        xlabel("% gait cycle. 0% and 100% = RHS")
    end

end
strides.normdata.Properties.VariableNames = varnames;


