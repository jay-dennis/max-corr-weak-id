function tabletotex(matin, rownames, colnames, outputname, Title, Caption)
%--------------------------------------------------------------------
% Jay Dennis, Feb 2016
%--------------------------------------------------------------------
% PURPOSE: output matlab data to formatted latex table (.tex)
%--------------------------------------------------------------------
% USAGE: tabletotex(matin, rownames, colnames, outputname, Title)
% where: matin    = data matrix (T x N)
%                   Note the 1st column must be the data type
%                   and the 2nd column must be the nlag
%      rownames   = cell array of row names,
%                   set = 0 if no row names
%      colnames   = cell array of column names in latex format
%      outputname = Output filename including directory path
%      Title      = Title (cell array) - each cell is a new line.
%      Caption    = Caption (string) (optional).
%--------------------------------------------------------------------
% RETURNS: 'outputname'.tex file into output directory.  Observe that you 
% must create the output directory if it does not exist.
% --------------------------------------------------------------------
% EXAMPLE:
%   clear table1 rownames colnames Title;
%   outputdir = 'output';
%   table1 = zeros(2,2);
%   rownames{1} = sprintf('Row 1 Name');
%   rownames{2} = sprintf('Row 2 Name');
%   colnames{1} = sprintf('Column 1 $(\\beta)$');
%   colnames{2} = sprintf('Col 2 $(\\sigma)$ ');
%   Title{1} = Title 1;
%   Title{2} = Title 2;
%   Caption = sprintf('Caption of Table');
%   outputname=sprintf('./%s/table_name', outputdir);  %Observe that the outputname 	
%						       %includes the output directory.
%   tabletotex(table1, rownames, colnames, outputname, Title);
% --------------------------------------------------------------------

fname = sprintf('%s.tex', outputname);

[T,N]=size(matin);

% Title = sprintf(' %s ', Title);
% Caption = sprintf(' %s ', Title);

%open file for writing
FID = fopen(fname, 'w');

% print the header info
fprintf(FID, ' \\begin{table}[H] \n');
fprintf(FID, ' \\small \n');
fprintf(FID, ' \\centering \n');
fprintf(FID, '\\begin{tabular}{|');
if iscell(rownames) == 0
    for j = 1:N
        fprintf(FID, 'c|');
    end;
    fprintf(FID, '} \n');
    for k = 1:length(Title)
        fprintf(FID, '\\multicolumn{%d}{c}{ %s } \\\\ \n', (N), Title{k});
    end
elseif iscell(rownames) ~= 0
    for j = 1:(N+1)
        fprintf(FID, 'c|');
    end;
    fprintf(FID, '} \n');
    for k = 1:length(Title)
        fprintf(FID, '\\multicolumn{%d}{c}{ %s } \\\\ \n', (N+1), Title{k});
    end
end;

%draw a horizontal line
fprintf(FID, ' \\hline \n');

%format the column string
colnames_string = '';
for a=1:length(colnames)
    colnames_string = [colnames_string '& ' colnames{a} ' '];
end
fprintf(FID, ' %s ', colnames_string);
fprintf(FID, ' \\\\ \n \\hline \n');

% print the table
for i = 1:T
  if iscell(rownames) ~= 0
    fprintf(FID, ' %s & ', rownames{i}); %Print Row Name in 1st Col
  end
  for j = 1:(N-1) %Print Table Row
    fprintf(FID, ' %4.4f & ', matin(i,j));
  end;
  fprintf(FID, ' %4.4f \\\\ \n', matin(i,N)); %Print last column
end;
fprintf(FID, ' \\hline \n');
fprintf(FID, '\\end{tabular}\n');

if nargin == 6
    fprintf(FID, ' \\caption{%s} \n', Caption);
end
fprintf(FID, ' \\end{table}\n');

% close the file and save
fclose(FID);
