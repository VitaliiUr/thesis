clear all

pkg load data-smoothing

filename="CROSS_old_30mev.dat";
column_number = 2;

data02=dlmread("./oldData/LO_30MeV.dat","",0,0);
data12=dlmread("./oldData/NLO_30MeV.dat","",0,0);
data22=dlmread("./oldData/N2LO_30MeV.dat","",0,0);
data32=dlmread("./oldData/N3LO_30MeV.dat","",0,0);
data42=dlmread("./oldData/N4LO_30MeV.dat","",0,0);


datar1=dlmread("./oldData/r08_30MeV.dat","",0,0);
datar2=dlmread("./oldData/r09_30MeV.dat","",0,0);
datar3=dlmread("./oldData/r1_30MeV.dat","",0,0);
datar4=dlmread("./oldData/r11_30MeV.dat","",0,0);
datar5=dlmread("./oldData/r12_30MeV.dat","",0,0);



data2(:,1)=data02(:,1);
data2(:,2)=data02(:,column_number);
data2(:,3)=data12(:,column_number);
data2(:,4)=data22(:,column_number);
data2(:,5)=data32(:,column_number);
data2(:,6)=data42(:,column_number);

data2(:,7)=datar1(:,column_number);
data2(:,8)=datar2(:,column_number);
data2(:,9)=datar3(:,column_number);
data2(:,10)=datar4(:,column_number);
data2(:,11)=datar5(:,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO\tNLO\tN2LO\tN3LO\tN4LO\tr0.8\tr0.9\tr1.0\tr1.1\tr1.2\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append', 'precision', 5);

for i = 1:size(data2)(1)
	deltaData(i,1) = data2(i,1);
	deltaData(i,2) = (max(data2(i,2:6))-min(data2(i,2:6)))/mean(data2(i,2:6));
	deltaData(i,3) = (max(data2(i,7:11))-min(data2(i,7:11)))/mean(data2(i,7:11));
endfor

deltaData(:,4) = regdatasmooth(deltaData(:,1),deltaData(:,2),"d",6);
deltaData(:,5) = regdatasmooth(deltaData(:,1),deltaData(:,3),"d",6);

deltaData(:,6) = abs(data2(:,6)-data2(:,5));
deltaData(:,7) = regdatasmooth(deltaData(:,1),deltaData(:,6),"d",3);

datfile = fopen("deltaData_old_30MeV.dat",'wt');
fprintf(datfile, 'DEG\tChiral\tCutoff\tChiral(smoothed)\tCutoff(smoothed)\tN3LO-N4LO\tN3LO-N4LO(smoothed)\n');
fclose(datfile)

dlmwrite("deltaData_old_30MeV.dat", deltaData,'delimiter', '\t', '-append', 'precision', 10);





clear all

filename="CROSS_old_100mev.dat";
column_number = 2;

data02=dlmread("./oldData/LO_100MeV.dat","",0,0);
data12=dlmread("./oldData/NLO_100MeV.dat","",0,0);
data22=dlmread("./oldData/N2LO_100MeV.dat","",0,0);
data32=dlmread("./oldData/N3LO_100MeV.dat","",0,0);
data42=dlmread("./oldData/N4LO_100MeV.dat","",0,0);


datar1=dlmread("./oldData/r08_100MeV.dat","",0,0);
datar2=dlmread("./oldData/r09_100MeV.dat","",0,0);
datar3=dlmread("./oldData/r1_100MeV.dat","",0,0);
datar4=dlmread("./oldData/r11_100MeV.dat","",0,0);
datar5=dlmread("./oldData/r12_100MeV.dat","",0,0);



data2(:,1)=data02(:,1);
data2(:,2)=data02(:,column_number);
data2(:,3)=data12(:,column_number);
data2(:,4)=data22(:,column_number);
data2(:,5)=data32(:,column_number);
data2(:,6)=data42(:,column_number);

data2(:,7)=datar1(:,column_number);
data2(:,8)=datar2(:,column_number);
data2(:,9)=datar3(:,column_number);
data2(:,10)=datar4(:,column_number);
data2(:,11)=datar5(:,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO\tNLO\tN2LO\tN3LO\tN4LO\tr0.8\tr0.9\tr1.0\tr1.1\tr1.2\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append', 'precision', 5);

for i = 1:size(data2)(1)
	deltaData(i,1) = data2(i,1);
	deltaData(i,2) = (max(data2(i,2:6))-min(data2(i,2:6)))/mean(data2(i,2:6));
	deltaData(i,3) = (max(data2(i,7:11))-min(data2(i,7:11)))/mean(data2(i,7:11));
endfor

deltaData(:,4) = regdatasmooth(deltaData(:,1),deltaData(:,2),"d",6);
deltaData(:,5) = regdatasmooth(deltaData(:,1),deltaData(:,3),"d",6);

deltaData(:,6) = abs(data2(:,6)-data2(:,5))
deltaData(:,7) = regdatasmooth(deltaData(:,1),deltaData(:,6),"d",3);


datfile = fopen("deltaData_old_100MeV.dat",'wt');
fprintf(datfile, 'DEG\tChiral\tCutoff\tChiral(smoothed)\tCutoff(smoothed)\tN3LO-N4LO\tN3LO-N4LO(smoothed)\n');
fclose(datfile)

dlmwrite("deltaData_old_100MeV.dat", deltaData,'delimiter', '\t', '-append', 'precision', 10);