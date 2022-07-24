clear all

filename="DeutAX_Deg_100mev.dat";
column_number = 5;

data02=dlmread("2.dat3.0.2.100p8mev","",1,0);
data12=dlmread("2.dat3.1.2.100p8mev","",1,0);
data22=dlmread("2.dat3.2.2.100p8mev","",1,0);
data32=dlmread("2.dat3.3.2.100p8mev","",1,0);
data41=dlmread("2.dat3.4.1.100p8mev","",1,0);
data42=dlmread("2.dat3.4.2.100p8mev","",1,0);
data43=dlmread("2.dat3.4.3.100p8mev","",1,0);
data44=dlmread("2.dat3.4.4.100p8mev","",1,0);
data52=dlmread("2.dat3.5.2.100p8mev","",1,0);
dataAV18=dlmread("2.dat3.av18.100mev","",1,0);
Pi=3.14159265359;

for i=1:length(data02)
  data02(i,1)=round(data02(i,1)*180/Pi);
%  data12(i,1)=round(data12(i,1)*180/Pi);
%  data22(i,1)=round(data22(i,1)*180/Pi);
%  data32(i,1)=round(data32(i,1)*180/Pi);
%  data41(i,1)=round(data41(i,1)*180/Pi);
%  data42(i,1)=round(data42(i,1)*180/Pi);
%  data43(i,1)=round(data43(i,1)*180/Pi);
%  data44(i,1)=round(data44(i,1)*180/Pi);
%  data52(i,1)=round(data52(i,1)*180/Pi);
%  dataAV18(i,1)=round(dataAV18(i,1)*180/Pi);
endfor

data2(:,1)=data02(:,1);
data2(:,2)=data02(:,column_number);
data2(:,3)=data12(:,column_number);
data2(:,4)=data22(:,column_number);
data2(:,5)=data32(:,column_number);
data2(:,6)=data41(:,column_number);
data2(:,7)=data42(:,column_number);
data2(:,8)=data43(:,column_number);
data2(:,9)=data44(:,column_number);
data2(:,10)=data52(:,column_number);
data2(:,11)=dataAV18(:,column_number);

datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO2\tNLO2\tN2LO2\tN3LO2\tN4LO1\tN4LO2\tN4LO3\tN4LO4\tN4LO+2\tAV18\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append');
