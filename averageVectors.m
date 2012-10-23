function [tVec,	M, A, Mmean, endVec, tfVec] = averageVectors(D,tD,rules)
% D 	- a matrix with rows of vectors. The combination of different length vectors into one matrix has to be done in the main program.
% tD 	- the corresponding matrix of time vectors for each data vector in D. there are no zeros in the time seris - so the zeros padding the end of the rows in the matrix is used to defer the size of the data set
% rules - a struct with variables used for averaging rules

% step 0
N = length(tD(:,1));
% step 1: find over lap of data
for i=1:N
	loc = find(tD(i,:)==0,1);
	if isempty(loc)
		endVec(i) = length(tD(1,:));
	else
		endVec(i) = find(tD(i,:)==0,1)-1;
        if endVec(i)<5 % arbitrary small number - to avoid problems with zero in the beginning of tD
            endVec(i) = find(tD(i,:)==0,1,'last')-1;
        end
	end
	tfVec(i) = tD(i,endVec(i));
    if tfVec(i)==0
        keyboard;
    end
end

tsVec = tD(:,1);
ts = max(tD(:,1));
tf = min(tfVec);

% step 3: interpolate the data to a common time vector
for i=1:N
	tStepVec(i) = min(diff(tD(i,1:endVec(i))));
end

if or(isempty(tStepVec),tStepVec==0)
    tStepVec = 1/24/60*10; % 10 minute default
end
tStepVec(find(tStepVec==0))=[];
tStep = min(tStepVec);
tVec = ts:tStep:tf;
for i=1:N
	M(i,:) = interp1(tD(i,1:endVec(i)),D(i,1:endVec(i)),tVec,'*nearest'); % the deafualt - "linear" method - makes an extra NaN for each "NaN hole" in the data
end

% step 2: find the holes in the data, and save a vector of the data with
% out holes - in every time step
% and normalize the data to <max = 1>

% well - the first part is easy. see nanmean.m.
if N==1
    Mmean=M;
    Mmean = Mmean / max(Mmean);
else
    % averaging
    Mmean = nanmean(M);
    
    % throwing outliers - because of problem with too many NaNs in the data the
    % averaging becomes nonsense if one data source has information and all the
    % others have NaN...
    % so outlier rule - is throwing data (that is - turning into NaN) if more
    % then half the stations have NaN
    Mmean(sum(isnan(M))>floor(N/2)) = NaN;
    
    % normalizing
    Mmean = Mmean / max(Mmean);
end
A = 1;

