close all;
clc;
clear all;
for i=1:4
    fprintf("value of i before loop= %.2f\n",i);
 if i==3
     continue;
 end
 fprintf("value of i AFTER loop= %.2f\n",i);
    fileid = ['ra',num2str(i),'.txt'];
 dlmread(fileid);
% a = [rafik, 3, 67ra];
%  i=1:length(a) 

%     a(i)
% i
end