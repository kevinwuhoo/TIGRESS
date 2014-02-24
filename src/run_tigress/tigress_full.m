function success=tigress_full(data,reg,varargin)

% Complete TIGRESS GRN inference method.
% Runs the TIGRESS gene regulatory network inference method given 
% expression data and a list of transcription factors. Options can be 
% included (see details below). It outputs a text file containing a 
% ranked list of edges with associated probabilities.
% For more information, please see :
% Haury et al.: 'TIGRESS: Trustful Inference of Gene REgulation using 
% Stability Selection',2012 
%
% Syntax 1: tigress_full(expression,tflist) % options to default
% Syntax 2: tigress_full(expression, tflist,'option',option_value,...)
% 
% REQUIRED INPUTS: 
%  - expression_file: a file containing the expression data for all genes. 
%  This input should be of the form 'path/expression_file_name'
%  - tflist_file: a file containing the list of transcription factors. This
%  input should be of the form 'path/expression_file_name'
%
% OPTIONAL INPUTS:
%  - R: number of resampling runs (default 1000)
%  - L: number of LARS steps at each iteration (default:5)
%  - alpha: randomization parameter alpha (in ]0,1]) (def:0.2)
%  - method: either 'original' or 'area' (default: 'area')
%  - cutoff: number of edges to predict (default: Inf)
%  - verbose: do you want the algorithm to show you the progress? 
%            (default=true)
%  - name_net: path/name.ext of file to write network into.
%              (Default='./edges.txt')
%
% OUTPUT:
% A file name 'edges.txt' that writes itself in the working directory.
% It contains 3 columns:
% - 1st column: transcription factors
% - 2nd column: target genes
% - 3rd column: probability of edge existence
% Note that the file is output such that the probabilities are ranked
% decreasingly.
% 
% Example:
% tigress_full(expression_file,tflist_file,'L',3,'R',500,'verbose',false)
%
% See also: tigress.m, score_edges.m, predict_network.m
%
% Anne-Claire Haury, 2012

success=-1;
%% Parse arguments
p = inputParser;   % Create an instance of the class.
p.addRequired('--data', @ischar);
p.addRequired('--reg', @ischar);
p.addParamValue('--L', '5', @ischar);
p.addParamValue('--R', '1000', @ischar);
p.addParamValue('--alpha', '.2', @ischar);
p.addParamValue('--method','area',@(x)any(strcmpi(x,{'area','original'})));
p.addParamValue('--cut', 'Inf', @ischar);
p.addParamValue('--verbose','1',@ischar);
p.addParamValue('--name_net',[data,'_TIGRESS_predictions.txt'],@isstr);
p.addParamValue('--LarsAlgo','lars',@(x)any(strcmpi(x,{'spams','glmnet','lars'})));
p.addParamValue('--parallel','0',@ischar);
p.parse('data','reg',varargin{:})



%% Show which arguments were not specified in the call.
disp(' ') 
disp 'List of arguments given default values:' 
for k=1:numel(p.UsingDefaults)
   field = char(p.UsingDefaults(k));
   value = num2str(p.Results.(field));
   if isempty(value)   
       value = '[]';   
   end
   fprintf('   ''%s''    defaults to %s \n', field, value)
end

%% Set variables
L=str2double(p.Results.L);
R=str2double(p.Results.R);
alpha=str2double(p.Results.alpha);
cut=str2double(p.Results.cut);
verbose=logical(str2double(p.Results.verbose));
parallel=logical(str2double(p.Results.parallel));
method=p.Results.method; 
LarsAlgo=p.Results.LarsAlgo;
name_net=p.Results.name_net;

%% Read inputs
data2=read_data(data,reg);

%% Get frequency matrix F
freq=tigress(data2,'L',L,'R',R,'alpha',alpha,'verbose',verbose,...
    'LarsAlgo',LarsAlgo,'parallel',parallel);

%% Get scores
scores=score_edges(freq,'method',method,'L',L);

%% Write edges
predict_network(scores,data2.tf_index,'genenames',data2.genenames,...
    'cut',cut,'name_net',name_net);
success=0;