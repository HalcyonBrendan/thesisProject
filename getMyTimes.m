function times=getMyTimes(Events,dirs)

%made by Greg, August 6 2011
% Extracts times of relevant strobes to timing structure from Sam's Events
% structure.

N=size(Events.strobes,1);
times.TrialStart=zeros(N,1);
times.T0=zeros(N,1);
times.T1=zeros(N,1);
times.T2=zeros(N,1);
times.T3=zeros(N,1);
times.T4=zeros(N,1);
times.T5=zeros(N,1);
times.T6=zeros(N,1);
times.T65=zeros(N,1);
times.T65sac=zeros(N,1); %monkey with these two s.t.: if not a saccade trial,
times.T66a=zeros(N,1);   %time of T65&T65sac are T6+600ms
times.T66on=zeros(N,1);
times.T66b=zeros(N,1);
times.T8=zeros(N,1);
times.Treach=zeros(N,1);
times.T9=zeros(N,1);
times.T10=zeros(N,1);
times.T13=zeros(N,1);
times.T18=zeros(N,1);
times.T19=zeros(N,1);


%% read strobes

for i=1:size(Events.strobes,1);
  try times.TrialStart(i)=Events.TrialStart(i);
  catch
    times.TrialStart(i)=Events.TrialStart0(i);
  end
  times.T0(i)=Events.strobeTimes(i,Events.strobes(i,1:10)==0)';
  times.T1(i)=Events.strobeTimes(i,Events.strobes(i,:)==1)';
  times.T2(i)=Events.strobeTimes(i,Events.strobes(i,:)==2)';
  times.T3(i)=Events.strobeTimes(i,Events.strobes(i,:)==3)';
  times.T4(i)=Events.strobeTimes(i,Events.strobes(i,:)==4)';
  times.T5(i)=Events.strobeTimes(i,Events.strobes(i,:)==5)';
  times.T6(i)=Events.strobeTimes(i,Events.strobes(i,:)==6)';
  times.T8(i)=Events.strobeTimes(i,Events.strobes(i,:)==8)';
  times.T9(i)=Events.strobeTimes(i,Events.strobes(i,:)==9)';
  times.T10(i)=Events.strobeTimes(i,Events.strobes(i,:)==10)';
  times.T13(i)=Events.strobeTimes(i,Events.strobes(i,:)==13)';
  
  %Times for 18 and 19 must be able to accomodate:
  %"18 19"
  %"518 18"
  %"50" "51"
  times.T18(i)=Events.strobeTimes(i,find(Events.strobes(i,:)==18 | Events.strobes(i,:)==518 | Events.strobes(i,:)==50,1) )';
  tmp=Events.strobeTimes(i,find(Events.strobes(i,:)==19 | Events.strobes(i,:)==519 | Events.strobes(i,:)==51,1))';
  if isempty(tmp)
    times.T19(i)=-1; %"518 18"
  else
    times.T19(i)=tmp;
  end
  if ~isempty(find(Events.strobes(i,:)==65, 1))
    times.T65(i)=Events.strobeTimes(i,Events.strobes(i,:)==65)';
  else
    times.T65(i)=nan; %watch out for this!
  end
  if ~isempty(find(Events.strobes(i,:)==66, 1))
    tmp=Events.strobeTimes(i,Events.strobes(i,:)==66)';
    times.T66a(i)=tmp(1)';
    times.T66b(i)=tmp(2)';
  else
    times.T66a(i)=nan;
    times.T66b(i)=nan;
  end
end

%%

%if provided with a directory for eye data, find time of 65 saccades
%if nargin>=2
  %try
    %load([dirs 'RecordingVariables.mat'])
    xeye=Eye.x;
    yeye=Eye.y;
    %hack solution to gain being high for this recording
    %if strcmp(dirs,'/Users/Mercy/Academics/MephistoArrays/MAD2/2010/12/06/M20101206_245/')
    %  xeye=xeye./100;
    %  yeye=yeye./100;
    %end
    for i=1:size(Events.strobes,1)
      if ~isempty(find(Events.strobes(i,:)==65, 1))
        tmp=myFindSaccade(xeye(i,:),yeye(i,:),times.T65(i),.07,4,Events.header(i,:))';
        times.T65sac(i)=tmp(1)';
      else
        times.T65sac(i)=nan; %watch out for this!
      end
    end
  %catch
  %end
%end

%% if provided with a directory for hand data, find time of reach
if nargin>=2
  try
    load([dirs 'RecordingVariables.mat'])
    xhand=Hand.x;
    yhand=Hand.y;
    for i=1:size(Events.strobes,1)
      tmp=myFindReach(xhand(i,:),yhand(i,:),times.T8(i))';
      if ~isempty(tmp)
        times.Treach(i)=tmp(1)';
      else
        times.Treach(i)=times.T8(i)';
      end
    end
  catch
  end
end

%% if provided with a directory for RecordingVariables, find time of pursuit onset

% two criteria for pursuit onset:  
%   position: z > 1
%   speed: speed > .25*avg_speed

%if nargin>=3
  %try
    %load([dirs 'RecordingVariables.mat'])
    
    b=fir1(12,.1); % FIR filter
    for m=1:size(Events.strobes,1)
      waitbar(m/size(Events.strobes,1))
      if ~isempty(find(Events.strobes(m,:)==66, 1))
        
        % get time
        t1=find(Eye.T>times.T6(m)+SPstart+50,1,'first');
        if abs(Events.header(m,6))>abs(Events.header(m,7)) % x or y direction?
          xx=Eye.x(m,t1-200:t1+1200)-mean(Eye.x(m,t1-200:t1));
        else
          xx=Eye.y(m,t1-200:t1+1200)-mean(Eye.y(m,t1-200:t1));
        end
        
        eyebl=xx(1:200);
        ff=filtfilt(b,1,double(xx));
        z=nan(801,1);
        eye_dir=sign(median(xx(end-50:end))-median(xx(1:50)));
        avgspeed=polyfit(1:1201,xx(201:end),1);
        avgspeed=abs(avgspeed(1));
        p=nan(801,2);
        for ti=51:801
          z(ti)=(ff(ti)-mean(eyebl))/std(eyebl)*eye_dir;
          p(ti,:)=polyfit(1:25,ff(ti-12:ti+12),1);
        end
        p(z>5,1)=10^6;
        tmp=strfind(z'>1 & abs(p(:,1))'>avgspeed*.25,ones(1,100));
        if ~isempty(tmp)
          T_puronset_rel66=tmp(1);
          times.T66on(m)=t1-200+T_puronset_rel66+Eye.T(1);
        else
          times.T66on(m)=nan;
        end
        
      else
        times.T66on(m)=nan; 
      end
    end
%  end
%end