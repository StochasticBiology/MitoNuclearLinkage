Rf=[0.05:0.05:0.95]
Rm=[0.05:0.05:0.95]
n=size(Rf,2);
Wf=zeros(1,n);
Wm=zeros(1,n);
dr=0.05;
for i=1:n
    i
    Wf(:,i)=NucMTPL(0.1,0.1,0.01,0.1,0.5, Rf(i),Rf(i)+dr, 0.05,1)
    %Wm(:,i)=NucMTPL(0.1,0.1,0.01,0.9,0.2, Rm(i),Rm(i)-dr, 0.05,1)
end
figure
hold
plot(Rf,Wf)
%plot(Rm,Wm)

% for i=1:n
%     i
%     Wf(:,i)=NucMTPL(0.1,0.1,0.01,0.9,0.1, Rf(i),Rf(i)+dr, 0.05,0)
%     %Wm(:,i)=NucMTPL(0.11,0.1,0.01,0.9,0.1, Rm(i),Rm(i)-dr, 0.05,0)
% end
% 
% plot(Rf,Wf)
% plot(Rm,Wm)

% for i=1:n
%     i
%     Wf(:,i)=NucMTPL(0.1,0.1,0.01,0.9,0.5, Rf(i),Rf(i)+dr, 0.05,0)
%     Wm(:,i)=NucMTPL(0.1,0.1,0.01,0.9,0.5, Rm(i),Rm(i)-dr, 0.05,0)
% end
% 
% plot(Rf,Wf)
% plot(Rm,Wm)