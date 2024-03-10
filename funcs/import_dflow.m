function [data, varnames] = import_dflow(dflow_file,dtype)

%% Imports data from an D-FLOW txt file into matlab
%
% INPUT VARIABLES:
%   dflow_file = full path for the marker or treadmill dflow file to be analyzed
%   dtype = "matrix" will output data all numeric matrix
%
% OUTPUT VARIABLES:
%   data = table (default). if dtype == "matrix", the data is all numeric data 
%   varnames = column headers, variable names
%
% written by Helen J. Huang (July 2023)
% updates tracked in git

file_id = fopen(dflow_file);

% Get column names
dfheader = fgetl(file_id);
varnames = string(regexp(dfheader, '\t', 'split'));
varnames = strrep(varnames,".","_");
N = length(varnames);

% to read in numeric data
formatspec = repmat("double",1,N);

fclose(file_id);

%% Read data in using 
dataLines = [2, Inf];

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", N);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = varnames;
opts.VariableTypes = formatspec;

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
datatable = readtable(dflow_file, opts);

if ~exist("dtype","var"), dtype = "default"; end 

if dtype == "matrix"
    data = table2array(datatable); 
else
    data = datatable;    
end

    