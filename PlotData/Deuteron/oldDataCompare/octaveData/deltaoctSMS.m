clear all

pkg load data-smoothing


% filename="CROSS_old_100mev.dat";
% column_number = 2;

% data02=dlmread("LO_100MeV.dat","",0,0);
% data12=dlmread("NLO_100MeV.dat","",0,0);
% data22=dlmread("N2LO_100MeV.dat","",0,0);
% data32=dlmread("N3LO_100MeV.dat","",0,0);
% data42=dlmread("N4LO_100MeV.dat","",0,0);


% datar1=dlmread("r08_100MeV.dat","",0,0);
% datar2=dlmread("r09_100MeV.dat","",0,0);
% datar3=dlmread("r1_100MeV.dat","",0,0);
% datar4=dlmread("r11_100MeV.dat","",0,0);
% datar5=dlmread("r12_100MeV.dat","",0,0);



% data2(:,1)=data02(:,1);
% data2(:,2)=data02(:,column_number);
% data2(:,3)=data12(:,column_number);
% data2(:,4)=data22(:,column_number);
% data2(:,5)=data32(:,column_number);
% data2(:,6)=data42(:,column_number);

% data2(:,7)=datar1(:,column_number);
% data2(:,8)=datar2(:,column_number);
% data2(:,9)=datar3(:,column_number);
% data2(:,10)=datar4(:,column_number);
% data2(:,11)=datar5(:,column_number);


% datfile = fopen(filename,'wt');
% fprintf(datfile, 'DEG\tLO\tNLO\tN2LO\tN3LO\tN4LO\tr0.8\tr0.9\tr1.0\tr1.1\tr1.2\n');
% fclose(datfile)

% dlmwrite(filename, data2,'delimiter', '\t', '-append', 'precision', 5);
data2 = dlmread("CROSS_30mev.dat","",1,0);
for i = 1:size(data2)(1)
	deltaData(i,1) = data2(i,1);
	deltaData(i,2) = (max([data2(i,2:5),data2(i,7),data2(i,10)])-min([data2(i,2:5),data2(i,7),data2(i,10)]))/mean([data2(i,2:5),data2(i,7),data2(i,10)]);
	deltaData(i,3) = (max(data2(i,6:9))-min(data2(i,6:9)))/mean(data2(i,6:9));
endfor

deltaData(:,4) = regdatasmooth(deltaData(:,1),deltaData(:,2),"d",5);
deltaData(:,5) = regdatasmooth(deltaData(:,1),deltaData(:,3),"d",5);

deltaData(:,6) = abs(data2(:,7)-data2(:,5));
deltaData(:,7) = regdatasmooth(deltaData(:,1),deltaData(:,6),"d",3);


datfile = fopen("deltaData_30MeV.dat",'wt');
fprintf(datfile, 'DEG\tChiral\tCutoff\tChiral(smoothed)\tCutoff(smoothed)\tN3LO-N4LO\tN3LO-N4LO(smoothed)\n');
fclose(datfile)

dlmwrite("deltaData_30MeV.dat", deltaData,'delimiter', '\t', '-append', 'precision', 10);



data2 = dlmread("CROSS_100mev.dat","",1,0);
for i = 1:size(data2)(1)
	deltaData(i,1) = data2(i,1);
	deltaData(i,2) = (max([data2(i,2:5),data2(i,7),data2(i,10)])-min([data2(i,2:5),data2(i,7),data2(i,10)]))/mean([data2(i,2:5),data2(i,7),data2(i,10)]);
	deltaData(i,3) = (max(data2(i,6:9))-min(data2(i,6:9)))/mean(data2(i,6:9));
endfor

deltaData(:,4) = regdatasmooth(deltaData(:,1),deltaData(:,2),"d",5);
deltaData(:,5) = regdatasmooth(deltaData(:,1),deltaData(:,3),"d",5);

deltaData(:,6) = abs(data2(:,7)-data2(:,5));
deltaData(:,7) = regdatasmooth(deltaData(:,1),deltaData(:,6),"d",3);


datfile = fopen("deltaData_100MeV.dat",'wt');
fprintf(datfile, 'DEG\tChiral\tCutoff\tChiral(smoothed)\tCutoff(smoothed)\tN3LO-N4LO\tN3LO-N4LO(smoothed)\n');
fclose(datfile)

dlmwrite("deltaData_100MeV.dat", deltaData,'delimiter', '\t', '-append', 'precision', 10);