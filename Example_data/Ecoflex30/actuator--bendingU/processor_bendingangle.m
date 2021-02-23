%w10_5c=[pixperonecm1066 VarName2];
%w10_5c0=[838 812];

t=w10_5c;
t0=w10_5c0;
theta=0;
Pressure=0;
time=[1/23:1/23:length(t)/23];

for i=1:1:length(t)
    theta(i)= 90+atand((t(i,2)-t0(2))/(t(i,1)-t0(1)));
    if (time(i)>=4 )
        Pressure(i)=(time(i)-4)/6+10;
    else
        Pressure(i)=0;
    end
end

hold on; plot(Pressure,theta,'k');
dlmwrite('8x8w10_5chamber.txt',[Pressure' theta']);
t=w7_5c;
t0=w7_5c0;
theta=0;
Pressure=0;
time=[1/23:1/23:length(t)/23];

for i=1:1:length(t)
    theta(i)= 90+atand((t(i,2)-t0(2))/(t(i,1)-t0(1)));
    if (time(i)>=4 )
        Pressure(i)=(time(i)-4)/6+10;
    else
        Pressure(i)=0;
    end
end

hold on; plot(Pressure,theta,'r');
dlmwrite('8x8w7_5chamber.txt',[Pressure' theta']);
t=w7_7c;
t0=w7_7c0;
theta=0;
Pressure=0;
time=[1/23:1/23:length(t)/23];

for i=1:1:length(t)
    theta(i)= 90+atand((t(i,2)-t0(2))/(t(i,1)-t0(1)));
    if (time(i)>=4 )
        Pressure(i)=(time(i)-4)/6+10;
    else
        Pressure(i)=0;
    end
end

hold on; plot(Pressure,theta,'g');
dlmwrite('8x8w7_7chamber.txt',[Pressure' theta']);
xlabel('Pressure(KPa)'); ylabel('Bending Angle (Degree)'); legend('8x4w10-5','8x4w7-5','8x4w7-7'); title('Input pressure Vs Angular displacement for Bendin SPAs');
