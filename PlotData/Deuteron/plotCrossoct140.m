clear all

filename="CROSS_140mev.dat";
column_number = 3;

data02=dlmread("data/2.dat1.0.2.140mev","",2,0);
data12=dlmread("data/2.dat1.1.2.140mev","",2,0);
data22=dlmread("data/2.dat1.2.2.140mev","",2,0);
data32=dlmread("data/2.dat1.3.2.140mev","",2,0);
data41=dlmread("data/2.dat1.4.1.140mev","",2,0);
data42=dlmread("data/2.dat1.4.2.140mev","",2,0);
data43=dlmread("data/2.dat1.4.3.140mev","",2,0);
data44=dlmread("data/2.dat1.4.4.140mev","",2,0);
data52=dlmread("data/2.dat1.5.2.140mev","",2,0);
%dataAV18=dlmread("data/2.dat1.av18.100mev","",2,0);


% Pi=3.14159265359;

for i=1:length(data02)
  data02(i,1)=round(data02(i,1)*180/pi);
endfor

data2(:,1)=data02(1:181,1);
data2(:,2)=data02(1:181,column_number);
data2(:,3)=data12(1:181,column_number);
data2(:,4)=data22(1:181,column_number);
data2(:,5)=data32(1:181,column_number);
data2(:,6)=data41(1:181,column_number);
data2(:,7)=data42(1:181,column_number);
data2(:,8)=data43(1:181,column_number);
data2(:,9)=data44(1:181,column_number);
data2(:,10)=data52(1:181,column_number);
%data2(:,11)=dataAV18(1:180,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO2\tNLO2\tN2LO2\tN3LO2\tN4LO1\tN4LO2\tN4LO3\tN4LO4\tN4LO+\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append');
