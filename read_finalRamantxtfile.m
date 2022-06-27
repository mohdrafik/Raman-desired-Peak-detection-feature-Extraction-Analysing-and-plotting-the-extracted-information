clc
close all;

fid = fopen('RamanData_file.txt','r');
%Advance five lines:
linesToSkip = 2;
for ii = 1:linesToSkip-1
    fgetl(fid);
end
%Process all remaining lines
tline =  fgetl(fid);
your_data = [22,8]; %You should allocate if you know how large your data is
while (~isempty(tline) )
    tline = fgetl(fid);
    %Getting rid of non-numbers
    tline = regexprep(tline,'[^0-9\s+-.eE]','');
    your_data = [your_data; str2num(tline)];
end
fclose(fid);