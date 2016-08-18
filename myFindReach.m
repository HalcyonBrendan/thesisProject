function reachtimes=myFindReach(eyeX,eyeY,t)

%made by Greg, March 9 2012.
% Calculates the time of the reach occuring after time 't', by finding a
% hand speed of >1 unit/ms.


reachtimes=nan;
if isnan(t)
  return
end
t=t+1000; %analog eye starts at -1000

if size(eyeX,1)==1 && size(eyeX,2)>1
  eyeX=eyeX';
end
if size(eyeY,1)==1 && size(eyeY,2)>1
  eyeY=eyeY';
end

speed=[diff(eyeX(t:end)) diff(eyeY(t:end))];
speed=sqrt(speed(:,1).^2+speed(:,2).^2); %combine x and y velocities
reachtimes=find(speed>1,1);
reachtimes=reachtimes+t*ones(size(reachtimes))-1000;
