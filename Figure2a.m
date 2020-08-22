R=[0.01:0.02:0.15,0.15:0.05:0.95,0.999]
n=size(R,2);
Wf=zeros(1,n);
Wm=zeros(1,n);
for i=1:n
    i
    Wf(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.9,0,0.5,0.1,R(i),0)
    Wm(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.9,0,0.5,0.1,R(i),1)
end
% figure
% hold
% box
plot(R,Wf)
plot(R,Wm)
% 
% for i=1:n
%     i
%     Wf(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.5,0,0.5,0,R(i),0)
%     Wm(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.5,0,0.5,0,R(i),1)
% end
% 
% plot(R,Wf)
% plot(R,Wm)
% 
% for i=1:n
%     i
%     Wf(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.1,0,0.5,0,R(i),0)
%     Wm(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.1,0,0.5,0,R(i),1)
% end
% 
% plot(R,Wf)
% plot(R,Wm)