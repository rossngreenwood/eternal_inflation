function test(outfile_name)

if nargin < 1
    outfile_name = ' eternal_inflation/Code/worker00.txt ';
end
outfile_name = outfile_name(2:end-1);

fid = fopen(outfile_name,'a');
if fid == -1
    disp('bad');
    fid = fopen(outfile_name,'w'); fclose(fid);
    fid = fopen(outfile_name,'a');
else
    disp('good');
end
fprintf(fid,'\r\nhello world');
fclose(fid);