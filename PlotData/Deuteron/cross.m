data=dlmread("data/CROSS_30mev.dat","",1,0);

data2(:,1)=data(:,1);
data2(:,2)=2*abs(data(:,2)-data(:,3))./(data(:,2)+data(:,3))*100;
data2(:,3)=2*abs(data(:,4)-data(:,3))./(data(:,4)+data(:,3))*100;
data2(:,4)=2*abs(data(:,5)-data(:,4))./(data(:,5)+data(:,4))*100;

datas=dlmread("data/CROSS_100mev.dat","",1,0);

datas2(:,1)=datas(:,1);
datas2(:,2)=2*abs(datas(:,4)-datas(:,3))./(datas(:,4)+datas(:,3))*100;
datas2(:,2)=2*abs(datas(:,2)-datas(:,3))./(datas(:,2)+datas(:,3))*100;
datas2(:,3)=2*abs(datas(:,4)-datas(:,3))./(datas(:,4)+datas(:,3))*100;
datas2(:,4)=2*abs(datas(:,5)-datas(:,4))./(datas(:,5)+datas(:,4))*100;
