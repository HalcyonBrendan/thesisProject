

% Make figure that shows percent truePosect, percent CCW1, OPP, CW1 for misses
% NOTE: 1,5=129=up, 2,6=130=down, 3,7=131=left, 4,8=132=right.

% Get some size constants
numIts = size(answerCheck,1);
maxTrainIdx = size(answerCheck,2);
numTests = size(answerCheck,3);

classes = unique(answerCheck(:,:,:,1));
numClasses = length(classes);

% Initialize some result counters
truePos = zeros(numClasses,1);
trueNeg = zeros(numClasses,1);
falsePos = zeros(numClasses,1);
falseNeg = zeros(numClasses,1);
ccw = zeros(numClasses,1);
cw = zeros(numClasses,1);
opp = zeros(numClasses,1);
tot = zeros(numClasses,1);
guessTot = zeros(numClasses,1);

for i = 1:numIts
    for j = 1:numTests
        % Find out what direction the look was towards
        if answerCheck(i,maxTrainIdx,j,1) == 1
            % Count total number that were supposed to be this direction
            tot(1) = tot(1)+1;
            % Find out what direction we guessed
            if answerCheck(i,maxTrainIdx,j,2) == 1
                truePos(1) = truePos(1)+1;
                guessTot(1) = guessTot(1)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 2
                opp(1) = opp(1)+1;
                guessTot(2) = guessTot(2)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 3
                ccw(1) = ccw(1)+1;
                guessTot(3) = guessTot(3)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 4
                cw(1) = cw(1)+1;
                guessTot(4) = guessTot(4)+1;
            end
        elseif answerCheck(i,maxTrainIdx,j,1) == 2
            tot(2) = tot(2)+1;
            if answerCheck(i,maxTrainIdx,j,2) == 1
                opp(2) = opp(2)+1;
                guessTot(1) = guessTot(1)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 2
                truePos(2) = truePos(2)+1;
                guessTot(2) = guessTot(2)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 3
                cw(2) = cw(2)+1;
                guessTot(3) = guessTot(3)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 4
                ccw(2) = ccw(2)+1;
                guessTot(4) = guessTot(4)+1;
            end            
        elseif answerCheck(i,maxTrainIdx,j,1) == 3
            tot(3) = tot(3)+1;
            if answerCheck(i,maxTrainIdx,j,2) == 1
                cw(3) = cw(3)+1;
                guessTot(1) = guessTot(1)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 2
                ccw(3) = ccw(3)+1;
                guessTot(2) = guessTot(2)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 3
                truePos(3) = truePos(3)+1;
                guessTot(3) = guessTot(3)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 4
                opp(3) = opp(3)+1;
                guessTot(4) = guessTot(4)+1;
            end            
        elseif answerCheck(i,maxTrainIdx,j,1) == 4
            tot(4) = tot(4)+1;
            if answerCheck(i,maxTrainIdx,j,2) == 1
                ccw(4) = ccw(4)+1;
                guessTot(1) = guessTot(1)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 2
                cw(4) = cw(4)+1;
                guessTot(2) = guessTot(2)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 3
                opp(4) = opp(4)+1;
                guessTot(3) = guessTot(3)+1;
            elseif answerCheck(i,maxTrainIdx,j,2) == 4
                truePos(4) = truePos(4)+1;
                guessTot(4) = guessTot(4)+1;
            end
        end
    end
end

for c = 1:numClasses
    falsePos(c) = guessTot(c)-truePos(c);
    falseNeg(c) = cw(c)+ccw(c)+opp(c);
    trueNeg(c) = numIts*numTests-truePos(c)-falsePos(c)-falseNeg(c);
end

rec(1) = truePos(1)/tot(1);
rec(2) = truePos(2)/tot(2);
rec(3) = truePos(3)/tot(3);
rec(4) = truePos(4)/tot(4);

prec(1) = truePos(1)/guessTot(1);
prec(2) = truePos(2)/guessTot(2);
prec(3) = truePos(3)/guessTot(3);
prec(4) = truePos(4)/guessTot(4);

recMic = sum(truePos)/sum(truePos+falseNeg);
recMac = sum(rec)/length(rec);
precMic = sum(truePos)/sum(truePos+falsePos);
precMac = sum(prec)/length(prec);

F1mic = 2*precMic*recMic/(precMic+recMic);
F1mac = 2*precMac*recMac/(precMac+recMac);


