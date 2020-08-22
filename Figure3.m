R=[0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.5]
n=size(R,2);
W=zeros(2,n);
for i=1:n
    i
    W(:,i)=NucMTRecomb(0.01,0.01,0.08,0.9,0.5, 0., 0.01, R(i))
end

plot(R,W(1,:))
plot(R,W(2,:))

% for i=1:n
%     i
%     W(:,i)=NucMTRecomb(0.1,0.1,0.01,0.9,0.5, 0., 0.01, R(i))
% end
% 
% plot(R,W(1,:))
% plot(R,W(2,:))
% 
% 
% for i=1:n
%     i
%     W(:,i)=NucMTRecomb(0.1,0.1,0.01,0.1,0.5, 0., 0.01, R(i))
% end
% 
% plot(R,W(1,:))
% plot(R,W(2,:))
