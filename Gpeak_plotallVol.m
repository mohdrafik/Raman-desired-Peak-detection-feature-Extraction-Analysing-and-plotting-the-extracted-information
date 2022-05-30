function y = Gpeak_plotallVol[datax,datay,i_file, n]
 
datay = (100*(i_file-1))*ones(size(datay))+ datay ;
figure
plot(datax,datay,'-ro',linewidth,'2',Marker,'b',MarkerSize,'3');
hold on 
end
