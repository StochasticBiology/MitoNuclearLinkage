Rf=[0.00:0.01:0.1, 0.1:0.1:0.48, 0.48];
%Rm=[0.02:0.05:0.1, 0.1:0.1:0.48, 0.5];
n=size(Rf,2);
Wf=zeros(1,n);
Wm=zeros(1,n);
dr=0.001;
for i=1:n
    i
    Wf(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.9,Rf(i),Rf(i)+dr,0,0.01,1)
    %Wm(:,i)=NucMTNewMutants(0.1,0.1,0.01,0.9,Rm(i),Rm(i)-dr,0,0.01,1)
end

hold
plot(Rf,Wf)
%plot(Rm,Wm)
% 
% for i=1:n
%     i
%     Wf(:,i)=NucMTNewMutants(0.12,0.10,0.01,0.1,R(i),R(i)+dr,0,0.1,1)
%     Wm(:,i)=NucMTNewMutants(0.12,0.10,0.01,0.1,R(i),R(i)-dr,0,0.1,1)
% end
% 
% plot(R,Wf)
% plot(R,Wm)