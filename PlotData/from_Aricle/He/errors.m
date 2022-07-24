
data0deg1=dlmread("./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-0deg","",9,0);
data0deg2=dlmread("./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-0deg","",9,0);
data0deg3=dlmread("./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-0deg","",9,0);
data0deg4=dlmread("./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-0deg","",9,0);

data60deg1=dlmread("./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-60deg","",9,0);
data60deg2=dlmread("./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-60deg","",9,0);
data60deg3=dlmread("./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-60deg","",9,0);
data60deg4=dlmread("./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-60deg","",9,0);

data120deg1=dlmread("./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-120deg","",9,0);
data120deg2=dlmread("./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-120deg","",9,0);
data120deg3=dlmread("./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-120deg","",9,0);
data120deg4=dlmread("./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-120deg","",9,0);

data180deg1=dlmread("./data/ppn-CUTNUM400-120.0-MeV_N4LO-incl-180deg","",9,0);
data180deg2=dlmread("./data/ppn-CUTNUM450-120.0-MeV_N4LO-incl-180deg","",9,0);
data180deg3=dlmread("./data/ppn-CUTNUM500-120.0-MeV_N4LO-incl-180deg","",9,0);
data180deg4=dlmread("./data/ppn-CUTNUM550-120.0-MeV_N4LO-incl-180deg","",9,0);


for i = 1:length(data0deg1)
	data0deg(i,1) = data0deg1(i,2);
	maxval = 10000*max([data0deg1(i,3),data0deg2(i,3),data0deg3(i,3),data0deg4(i,3)]);
	minval = 10000*min([data0deg1(i,3),data0deg2(i,3),data0deg3(i,3),data0deg4(i,3)]);
	avval = 10000*mean([data0deg1(i,3),data0deg2(i,3),data0deg3(i,3),data0deg4(i,3)]);
	if (avval > 1e-3)
		data0deg(i,2) = (maxval-minval)/avval;
	else
		data0deg(i,2) = 0;
	endif
endfor


for i = 1:length(data60deg1)
	data60deg(i,1) = data60deg1(i,2);
	maxval = 10000*max([data60deg1(i,3),data60deg2(i,3),data60deg3(i,3),data60deg4(i,3)]);
	minval = 10000*min([data60deg1(i,3),data60deg2(i,3),data60deg3(i,3),data60deg4(i,3)]);
	avval = 10000*mean([data60deg1(i,3),data60deg2(i,3),data60deg3(i,3),data60deg4(i,3)]);
	if (avval > 1e-3)
		data60deg(i,2) = (maxval-minval)/avval;
	else
		data60deg(i,2) = 0;
	endif
endfor


for i = 1:length(data120deg1)
	data120deg(i,1) = data120deg1(i,2);
	maxval = 10000*max([data120deg1(i,3),data120deg2(i,3),data120deg3(i,3),data120deg4(i,3)]);
	minval = 10000*min([data120deg1(i,3),data120deg2(i,3),data120deg3(i,3),data120deg4(i,3)]);
	avval = 10000*mean([data120deg1(i,3),data120deg2(i,3),data120deg3(i,3),data120deg4(i,3)]);
	if (avval > 1e-3)
		data120deg(i,2) = (maxval-minval)/avval;
	else
		data120deg(i,2) = 0;
	endif
endfor


for i = 1:length(data180deg1)
	data180deg(i,1) = data180deg1(i,2);
	maxval = 10000*max([data180deg1(i,3),data180deg2(i,3),data180deg3(i,3),data180deg4(i,3)]);
	minval = 10000*min([data180deg1(i,3),data180deg2(i,3),data180deg3(i,3),data180deg4(i,3)]);
	avval = 10000*mean([data180deg1(i,3),data180deg2(i,3),data180deg3(i,3),data180deg4(i,3)]);
	if (avval > 1e-3)
		data180deg(i,2) = (maxval-minval)/avval;
	else
		data180deg(i,2) = 0;
	endif
endfor

figure(1)
plot(data0deg(:,1),data0deg(:,2),"*")
title("0deg")
figure(2)
plot(data60deg(:,1),data60deg(:,2),"*")
title("60deg")
figure(3)
plot(data120deg(:,1),data120deg(:,2),"*")
title("120deg")
figure(4)
plot(data180deg(:,1),data180deg(:,2),"*")
title("180deg")
