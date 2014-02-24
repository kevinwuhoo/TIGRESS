function data=read_data(networkname)

% Transforms data from files into a structure
% 
% Syntax: data=read_data(expression_file,tflist_file)
% 
% This program assumes that: 
%  - your data is stored in ../data/<networkname>_Data
%  - all files within this folder begin with <networkname>
%  - this folder contains at least '<networkname>_expression_data.tsv'
%  - optionally, it may also contain '<neworkname>_transcription_factors.tsv', 
%  '<networkname>_gold_standard.tsv' and AUPR/AUROC distributions (see, e.g. DREAM5 data)
%  
%
% REQUIRED INPUTS: 
%  - networkname: string containing the name of your network
%
% See also: tigress_full
%
% Anne-Claire Haury, 2012

fprintf('Looking for data related to %s ... \n', networkname)

%% Get all files
datapath=['../data/',networkname,'_Data/'];
expression_file=[datapath,networkname,'_expression_data.tsv'];

fprintf('Getting expression data ... \n')
try 
      A=importdata(expression_file); 

catch err
      if ~exist(datapath,'dir')
          currentpath=pwd;
          [pathstr currentdir]=fileparts(currentpath);
          if any(strcmp(currentdir,{'data','doc','src'}))
              msg=sprintf('\n Cannot find directory. Try another network name? ');
              reply=input(msg,'s');
              if any(strcmp(reply,{'y',''}))
                  newname= input('\n Enter new name : ','s');
              else
                  throw(err);
              end
              expression_file2=['../data/',newname,'_Data/',newname,'_expression_data.tsv'];
              try 
                  A=importdata(expression_file2);
              catch err2
                  err3=addCause(err2,err);
                  throw(err3)
              end
              datapath=['../data/',newname,'_Data/'];
              networkname=newname;
          end
          
      end
end
genenames=A.textdata;
expdata=A.data;

%% Set data structure
data.expdata=expdata;
data.genenames=genenames';


%% Add optional fields
GSname=[datapath,networkname,'_gold_standard.tsv'];
pdfAUROCname=[datapath,networkname,'_AUROC.mat'];
pdfAUPRname=[datapath,networkname,'_AUPR.mat'];
tflist_file=[datapath,networkname,'_transcription_factors.tsv'];
if exist(tflist_file,'file')
   fprintf('Getting transcription factors ... \n')
   tflist=importdata(tflist_file);
   data.tflist=tflist;
   data.tf_index=find(ismember(data.tflist,data.genenames));
else
   fprintf('Generating transcription factors : all genes ... \n')
   data.tflist=genenames;
   data.tf_index=1:length(genenames);
end
if exist(GSname,'file')
    fprintf('Getting gold standard data ... \n')
    B = importdata(GSname);     
    idx=logical(B.data);
    tmp=B.data(idx);
    GS=B.textdata(idx,:);
    GSnum=zeros(length(tmp),3);
    GSnum(:,3)=tmp;
    N = size(GS,1);
    C = char(GS);
    D = C(:,2:end);
    F = str2num(D);
    GSnum(:,1:2) = [ F(1:N) F(N+1:2*N) ];
    data.gold_edges=GS;
    data.gold_edges_num=GSnum;
end
if exist(pdfAUROCname,'file')
    fprintf('Getting AUROC data... \n')
    pdf_auroc=load(pdfAUROCname);
    data.pdf_auroc=pdf_auroc;
end
if exist(pdfAUPRname,'file')
    fprintf('Getting AUPR data ... \n')
    pdf_aupr=load(pdfAUPRname);
    data.pdf_aupr=pdf_aupr;
end