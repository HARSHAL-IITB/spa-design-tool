% Processing of the 8x8w9 various samples
%Need to load following data to workspace
%Old_8x8w9;
%New8x8w9s1;
clc;
w7_5c=dlmread('8x8w7_5chambers.txt');
w7_7c=dlmread('8x8w7_7chambers.txt');
w10_5c=dlmread('8x8w10_5chambers.txt');
t=w7_5c;%Old_8x8w9;
Avg_ground=0;
Avg_count1=0;
Avg_dc=0;
Avg_count2=0;
state_change=0;
prev_state=0;
current_state=0;
Avg=zeros(1,3);
j=1;
for(i=1:1:length(t))
    if (t(i,1)>3 & t(i,3)==0)
        Avg_ground=Avg_ground + t(i,2);
        Avg_count1=Avg_count1+1;


    elseif(1.5<t(i,1)<3 & t(i,3)~=0)
        Avg_dc=Avg_dc+t(i,2);
        Avg_count2=Avg_count2+1;

    end
    current_state=t(i,3);

    if(prev_state~=current_state)
        if(current_state==0)
            Avg(j,2) =Avg_dc/Avg_count2;
            Avg(j,1) =prev_state;
            j=j+1;
        else
            Avg(j,3) =Avg_ground/Avg_count1;
        end
        Avg_ground=0;
        Avg_count1=0;
        Avg_dc=0;
        Avg_count2=0;   

    end

    prev_state=current_state;
end
Avg=sortrows(Avg)
Avg_w7_5c=[Avg(:,1) Avg(:,2)]

t=w7_7c;%New8x8w9s1;
Avg_ground=0;
Avg_count1=0;
Avg_dc=0;
Avg_count2=0;
state_change=0;
prev_state=0;

current_state=0;
Avg=zeros(1,3);
j=1;
for(i=1:1:length(t))
    if (t(i,1)>3 & t(i,3)==0)
        Avg_ground=Avg_ground + t(i,2);
        Avg_count1=Avg_count1+1;


    elseif(1.5<t(i,1)<3 & t(i,3)~=0)
        Avg_dc=Avg_dc+t(i,2);
        Avg_count2=Avg_count2+1;

    end
    current_state=t(i,3);

    if(prev_state~=current_state)
        if(current_state==0)
            Avg(j,2) =Avg_dc/Avg_count2;
            Avg(j,1) =prev_state;
            j=j+1;
        else
            Avg(j,3) =Avg_ground/Avg_count1;
        end
        Avg_ground=0;
        Avg_count1=0;
        Avg_dc=0;
        Avg_count2=0;   

    end

    prev_state=current_state;
end
Avg=sortrows(Avg)
Avg_w7_7c=[Avg(:,1) Avg(:,2)]

t=w10_5c;%New8x8w9s1;
Avg_ground=0;
Avg_count1=0;
Avg_dc=0;
Avg_count2=0;
Std_dc=0;
Std_ground=0;
state_change=0;
prev_state=0;

current_state=0;
Avg=zeros(1,4);
j=1;
for(i=1:1:length(t))
    if (t(i,1)>3 & t(i,3)==0)
        Avg_ground=Avg_ground + t(i,2);
        Avg_count1=Avg_count1+1;
        Std_ground=[Std_ground; t(i,2)];

    elseif(1.5<t(i,1)<3 & t(i,3)~=0)
        Avg_dc=Avg_dc+t(i,2);
        Avg_count2=Avg_count2+1;
        Std_dc=[Std_dc; t(i,2)]; 
    end
    current_state=t(i,3);

    if(prev_state~=current_state)
        if(current_state==0)
            Avg(j,2) =Avg_dc/Avg_count2;
            Avg(j,1) =prev_state;
            Std_dc
            Avg(j,4) =std(Std_dc);
            j=j+1;
        else
            Avg(j,3) =Avg_ground/Avg_count1;
        end
        Avg_ground=0;
        Avg_count1=0;
        Avg_dc=0;
        Avg_count2=0;   
        Std_dc=0;
        Std_ground=0;
    end
    prev_state=current_state;
end
Avg=sortrows(Avg)
Avg_w10_5c=[Avg(:,1) Avg(:,2) Avg(:,4)]



figure(2);
%hold on; plot(Avg_w7_4c(:,1),Avg_w7_4c(:,2),'*-')
hold on; plot(Avg_w7_5c(:,1),Avg_w7_5c(:,2),'*-g'); 
hold on; plot(Avg_w7_7c(:,1),Avg_w7_7c(:,2),'*-r')
hold on; plot(Avg_w10_5c(:,1),Avg_w10_5c(:,2),'*-k');errorbar(Avg_w10_5c(:,2),Avg_w10_5c(:,3));
hold off;
xlabel('Pressure(N/m^2)');
ylabel('Bending Force(N)');
title('Bending force for bending SPA 8x8w7 and 8x8w10 comparision');
legend ('w7_5chamber','w7_7chamber','w10_5chamber');
