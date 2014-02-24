function data=read_data(expression_file,tflist_file)

% Transforms data from files into a structure
% 
% Syntax: data=read_data(expression_file,tflist_file)
% 
% REQUIRED INPUTS: 
%  - expression_file: a file containing the expression data for all genes. 
%  This input should be of the form 'path/expression_file_name'
%  - tflist_file: a file containing the list of transcription factors. This
%  input should be of the form 'path/expression_file_name'
%
% See also: tigress_full
%
% Anne-Claire Haury, 2012

%% Import expression data
A=importdata(expression_file);
genenames=A.textdata;
expdata=A.data;
tflist=importdata(tflist_file);

%% Set data structure
data.expdata=expdata;
data.genenames=genenames;
data.tflist=tflist;
data.tf_index=find(ismember(data.tflist,data.genenames));