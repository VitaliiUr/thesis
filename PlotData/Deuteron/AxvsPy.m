dataAx=dlmread("data/DeutAX_Deg_100mev.dat","",1,0);
dataPy=dlmread("data/ProtonPolDeg_100mev.dat","",1,0);
% data30r = data30(:,1);
for i = 1:180
	% if (abs(mean([dataAx(i,7),dataAx(i,10)])) > 1e-2)
	% 	dataAxChiral(i) = abs(dataAx(i,7)-dataAx(i,10))/mean([dataAx(i,7),dataAx(i,10)]);
	% else
	% 	dataAxChiral(i) = 0;
	% endif
	if (abs(mean(dataAx(i,6:9))) > 1e-2)
		dataAxCut(i) = (max(dataAx(i,6:9))-min(dataAx(i,6:9)))/abs(mean(dataAx(i,6:9)));
	else
		dataAxCut(i) = 0;
	endif
	% if (abs(mean([dataPy(i,7),dataPy(i,10)])) > 2e-2)
	% 	dataPyChiral(i) = abs(dataPy(i,10)-dataPy(i,7))/mean([dataPy(i,7),dataPy(i,10)]);
	% else
	% 	dataPyChiral(i) = 0;
	% endif
	if (abs(mean(dataPy(i,6:9))) > 1e-2)
		dataPyCut(i) = (max(dataPy(i,6:9))-min(dataPy(i,6:9)))/abs(mean(dataPy(i,6:9)));
	else
		dataPyCut(i) = 0;
	endif
endfor

[el,num]=max(dataPy(:,10));
Pydif = abs(dataPy(num,10)-dataPy(num,7))/abs(mean([dataPy(num,10),dataPy(num,7)]));
cstrcat("difference for Py in point of max, theta = ",num2str(dataPy(num,1)),", = ",num2str(Pydif*100)," %")
cstrcat("difference for Py in point of max wrt cut-off, theta = ",num2str(dataPy(num,1)),", = ",num2str(dataPyCut(num)*100)," %")

[el,num]=min(dataPy(:,10));
Pydif = abs(dataPy(num,10)-dataPy(num,7))/abs(mean([dataPy(num,10),dataPy(num,7)]));
cstrcat("difference for Py in point of min, theta = ",num2str(dataPy(num,1)),", = ",num2str(Pydif*100)," %")
cstrcat("difference for Py in point of min wrt cut-off, theta = ",num2str(dataPy(num,1)),", = ",num2str(dataPyCut(num)*100)," %")

[el,num]=min(dataAx(:,10));
Axdif = abs(dataAx(num,10)-dataAx(num,7))/abs(mean([dataAx(num,10),dataAx(num,7)]));
cstrcat("difference for Ax in point of min, theta = ",num2str(dataAx(num,1)),", = ",num2str(Axdif*100)," %")
cstrcat("difference for Ax in point of min wrt cut-off, theta = ",num2str(dataAx(num,1)),", = ",num2str(dataAxCut(num)*100)," %")
