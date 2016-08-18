function saccadetimes=myFindSaccade(eyeX,eyeY,t,thresh,algflag,header)

%made by Greg, November 7 2011.
% Calculates the start and finish of the major saccades occuring after
% time 't'. Eye data and t are in ms, t counts from the beginning
% of eye data (first data point is 1 ms).
% algflag=1 - findSequential algorithm (legacy)
% algflag=2 - findJump algorithm


if nargin<4
  thresh=.3;
end
if nargin<5
  algflag=2;
end

saccadetimes=nan;
if isnan(t)
  return
end
t=t+1001; %analog eye starts at -1000

if size(eyeX,1)==1 && size(eyeX,2)>1
  eyeX=eyeX';
end
if size(eyeY,1)==1 && size(eyeY,2)>1
  eyeY=eyeY';
end

if algflag==1 % findSequential
  speed=[smooth(diff(eyeX(t:end)),5) smooth(diff(eyeY(t:end)),5)];
  speed=sqrt(speed(:,1).^2+speed(:,2).^2); %combine x and y velocities
  saccadetimes=findSequential(speed>thresh,10,2);
  while isempty(saccadetimes)
    thresh=thresh*.75;
    saccadetimes=findSequential(speed>thresh,10,2);
  end
  % remove saccades within twait=100 of each other
  dind=diff(saccadetimes);
  saccadetimes=saccadetimes([1 find(dind>100)+1]);
  saccadetimes=saccadetimes+t-1001;
  
elseif algflag==2 % findJump
  
  speed=[smooth(diff(smooth(eyeX(t:end)))) smooth(diff(smooth(eyeY(t:end))))];
  speed=sqrt(speed(:,1).^2+speed(:,2).^2); %combine x and y velocities
  eye=smooth(sqrt(eyeX(t:end).^2+eyeY(t:end).^2));

  tmp2=zeros(size(speed));
  for r=10:length(speed)
    tmp2(r)=sum(speed(r-9:r));
  end
  
  %[~,I]=sort(tmp2);
  %thresh=std(tmp2(I(1:round(length(tmp2)*.99))));
  thresh=std(tmp2(tmp2<1));
  
  sactimes=find(tmp2>thresh)-9;
  
  % remove sac times w/in 30ms of each other
  tmp=diff([0;sactimes]); 
  sactimes=sactimes(tmp>30);
  
  % remove sac times in the first and last 50 ms for convenience 
  % (these don't matter anyway)
  sactimes=sactimes(sactimes>50 & sactimes<length(speed)-50);
  
  % remove sac times where the eye doesn't change by at least 0.15
  tmp=abs(eye(sactimes+30)-eye(sactimes-10));
  sactimes=sactimes(tmp>.07);
  
  saccadetimes=sactimes+t-1001;
  
elseif algflag==3
  
  saccadetimes=findJump(eye,20,thresh);
  saccadetimes=saccadetimes+t-1001;
  
elseif algflag==4
  
  b=fir1(12,.1); % FIR filter
  if abs(header(6))>abs(header(7)) % x or y direction?
    xx=eyeX-mean(eyeX(t-200:t));
  else
    xx=eyeY-mean(eyeY(t-200:t));
  end
  
  ff=filtfilt(b,1,double(xx))';
  ff_snip=ff(t-200:t+1000); % take 200 ms before T65 to 1000 ms after
  avgspeed=abs(polyfit(1:length(ff_snip),ff_snip,1)); % get the speed over the whole bl
  speed=smooth(diff(ff_snip),2); % speed
  tmp=strfind(abs(speed)'>avgspeed(1)*5,ones(1,20));
  tmp=tmp(1)+5;
  saccadetimes=tmp+t-1001-200;

else 
  error('No algorithm selected')
end



function ind=findSequential(inp,N,wiggle)
%In input 'inp' find N sequential points (i.e. separated by 1 index),
%allowing for 'wiggle' points to be out. 'inp' is integer indices. Returns
%the first index in the found series.
%For example, inp=[1 2 3 4 6 7 9], N=9, wiggle=2 -> ind=1.
%Finds all series in 'inp' separated by 100 indices. This post-saccadic
%refractory period might need refinement.
% New algorithm: swipe 9501x1 with a filter, wait for product>N-wiggle
c=1;ind=[];twait=100;
for i=1:length(inp)
  r=min(i,N);
  filter=ones(r,1);
  pr=inp(i-r+1:i).*filter;
  if sum(pr)>N-wiggle
    ind(c)=i-round(r/2);
    c=c+1;
  end
end

% for i=1:length(inp)
%   r=min(i,N+1)-1;
%   % is there a string of 
%   if length(find(diff(inp(i-r:i))==1))>=N-wiggle-1 && ~twait
%     ind(c)=inp(i-round(r/2));
%     c=c+1;
%     %now wait 100ms until you can start the next saccade
%     twait=inp(i)+100;
%   end
%   if twait-inp(i)<0
%     twait=0;
%   end
% end



function ind=findJump(inp,N,jump)
%New algorithm. Finds saccades by looking at two points N ms apart. If the
%difference in y-value between points is >jump, a saccade occurred.
% ie. if dy=y(i)-y(i-N)>jump, keep i
%Control. If dy2=mean(y(i:i+20))-y(i-N)<jump, throw away i

inp2=inp(N+1:end);
inp1=inp(1:end-N);
dinp=abs(inp2-inp1);
ind=find(dinp>jump);
dind=diff(ind); % remove sac times w/in 100ms of each other
ind=ind(dind>100);
