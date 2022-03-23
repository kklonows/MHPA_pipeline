%% Section refnt
more on
format shortG

% import table in mergez.txt file as a variable called mergez
chrNum = dlmread('chrnum.txt');

if(chrNum == 23);
    mergez = dlmread('newMergez.txt');
else
    mergez = dlmread('mergez.txt');
end;

% import the variable values that were stored in the
% matlab_input_genomeADcoord.mat file

load('matlab_input_genomeADcoord.mat');

%need to generate refnt file for TSC1, TSC2, done manually.

% Load the minimum allele frequency that was inputted by the user in
% the python part of the pipeline
minfreq = dlmread('minfreq.txt');

% Load the minimum indel frequency that was inputted by the user in the
% python part of the pipeline
minIndelFreq = dlmread('minIndelFreq.txt');

% Load the minimum read count that was inputted by the user in the python
% part of the pipeline

minreadcount = dlmread('minreadcount.txt');

%{
refnt=zeros(137002,2);
refnt(1:75001,1)=9;
refnt(75002:137002,1)=16;
for i=1:75001;
    refnt(i,2)=135756000+i-1;
end;
for i=75002:137002;
    refnt(i,2)=2087000+i-75002;
end;
%}

initIndex = dlmread('initindex.txt');
endIndex = dlmread('endindex.txt');

% TXNIP chr1:145438462-145442628
seqLength = endIndex - initIndex + 1;
refnt=zeros(seqLength, 2);
refnt(1:seqLength,1)=chrNum;
for i=1:seqLength;
    refnt(i,2)=initIndex+i-1;
end;

exon_coord_flush = dlmread('exonCoordData.txt');

TSC1_exon_coord_flush = exon_coord_flush;
TSC2_exon_coord_flush = exon_coord_flush;
T1T2exonsflush = exon_coord_flush;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %FLCN chr17:17113000-17143000
% refnt=zeros(30000,2);
% refnt(1:30000,1)=17;
% for i=1:30000;
%     refnt(i,2)=17113000+i-1;
% end;

% %PIK3CA
% refnt=zeros(91979,2);
% refnt(1:91979,1)=3;
% for i=1:91979;
%     refnt(i,2)=178865902+i-1;
% end;
%
% %MTOR
% refnt=zeros(155972,2);
% refnt(1:155972,1)=1;
% for i=1:155972;
%     refnt(i,2)=11166592+i-1;
% end;
%
% %AKT1
% refnt=zeros(26402,2);
% refnt(1:26402,1)=14;
% for i=1:26402;
%     refnt(i,2)=105235686+i-1;
% end;
%
% %AKT2
% refnt=zeros(55219,2);
% refnt(1:55219,1)=19;
% for i=1:55219;
%     refnt(i,2)=40736224+i-1;
% end;
%
% % AKT3
% refnt=zeros(362846,2);
% refnt(1:362846,1)=1;
% for i=1:362846;
%     refnt(i,2)=243651535+i-1;
% end;
%
% % PTEN
% refnt=zeros(108817,2);
% refnt(1:108817,1)=10;
% for i=1:108817;
%     refnt(i,2)=89622870+i-1;
% end;

% %DEPDC5
% refnt=zeros(153083,2);
% refnt(1:153083,1)=22;
% for i=1:153083;
%     refnt(i,2)=32149937+i-1;
% end;
%% Section mergez-mergeza
% Set the number of samples to be the number of bam files in the directory
s = dlmread('numSamples.txt');

% s = the number of sequence files being analyzed

% Optimization: When the mergez array is traversed from 1 to s, then the
% variable for the array appears to change size on every loop iteration. To
% avoid this, I changed the order of traversal of the mergez array so that
% the array grows along the last dimension. So, the first time the mergez
% array is assigned a value, it is assigned to the element at position s
% and thus, Matlab allocates space for all of the s elements of the mergez
% array.
for x=s:-1:1;
    mergez(:,(x-1)*13+12)=0;
    mergez(:,(x-1)*13+13)=0;
end;
%%
%%
% mergez(:,~any(mergez,1))=[];w
% %need to modify merged call data from python output -  using pico - simply remove first line (of text)
% !!!need to import merged file, called all.b, and add in a line of zeros
%%

%below aligns all the chr nt to the same line, since are not aligned in input
%mergez=all;
%note that here we are assuming that input - mergez - had 13 columns per
%chr:nt:sample
mergeza=zeros(size(refnt,1),s*13);
for h=1:s;
    j=1;
    for i=1:size(refnt,1);
        if refnt(i,1)==mergez(j,(h-1)*13+1) && refnt(i,2)==mergez(j,(h-1)*13+2);
            for k=1:13;
                mergeza(i,(h-1)*13+k)=mergez(j,(h-1)*13+k);
            end;
            j=j+1;
        else
        end;
    end;
end;

% %rows in mergeza sequentially go through all nt positions in target genes
% %each sample gets 11 columns; entries are present only if there is some read data, otherwise chr and nt are 0, as are all other entries.  11 columns are: chr nt Af Gf Cf Tf Ar Gr Cr Tr indels
% %check alignment in a few rows:
% mergeza(10001:10002,1:451)
% mergeza(100001:100002,1:451)
% mergeza(1:2,1:451)

%% Section worka-workaa -> could be done earlier
%now delete empty rows:
%and eliminate 2 extra columns (12 and 13) for each sample:
worka=zeros(size(refnt,1),s*11);
k=0;
for i=1:size(refnt,1);
    if max(mergeza(i,1:s*13))>0;
        k=k+1;
        for j=1:s;
            for m=1:11;
            worka(k,(j-1)*11+m)=mergeza(i,(j-1)*13+m);
            end;
        end;
    else
    end;
end;
%k;

%k=78466; now delete rows with < 2 reads per sample
%rowinf used here to track which sample had how many reads, nt position.
workaa=zeros(k,s*11);
rowinf=zeros(100000,3);
h=0;
for i=1:k;
    sum=0;
    for j=1:s;
        if (worka(i,(j-1)*11+3)+worka(i,(j-1)*11+4)+worka(i,(j-1)*11+5)+worka(i,(j-1)*11+6)+worka(i,(j-1)*11+7)+worka(i,(j-1)*11+8)+worka(i,(j-1)*11+9)+worka(i,(j-1)*11+10)+worka(i,(j-1)*11+11))>1
            sum=1;
            maxreads=worka(i,(j-1)*11+3)+worka(i,(j-1)*11+4)+worka(i,(j-1)*11+5)+worka(i,(j-1)*11+6)+worka(i,(j-1)*11+7)+worka(i,(j-1)*11+8)+worka(i,(j-1)*11+9)+worka(i,(j-1)*11+10)+worka(i,(j-1)*11+11);
            sample=j;
        else
        end;
    end;
    if sum>0;
        h=h+1;
        for j=1:s*11;
            workaa(h,j)=worka(i,j);
        end;
        rowinf(h,1)=sample;
        rowinf(h,2)=maxreads;
        rowinf(h,3)=workaa(h,2);
    else
    end;
end;
% h;
%h=58492 at this point
%rows in workaa contain all nt positions in TSC1, then TSC2, with some reads
%each sample gets 11 columns; entries are present only if there is some read data, otherwise chr and nt are 0, as are all other %entries.  11 columns are: chr nt Af Gf Cf Tf Ar Gr Cr Tr indels

%%
%put chr and nt in 1st positions
for i=1:h;
    for j=1:s-1;
        workaa(i,1)=max(workaa(i,1),workaa(i,j*11+1));
        workaa(i,2)=max(workaa(i,2),workaa(i,j*11+2));
    end;
end;

mergezb=workaa;

%%
%now capture exon variants only +-25 nt
% mergezb=zeros(h,s*11);
% hh=h;
% h=1;
% for i=1:hh;
% for j=1:62;
% if ((workaa(i,1)==exons(j,1))&&(workaa(i,2)>(exons(j,2)-25))&&(workaa(i,2)<(exons(j,3)+25)));
% 	for k=1:s*11;
%        mergezb(h,k)=workaa(i,k);
%      end;
%      h=h+1;
% else
%   end;
% end;
% end;
% h=h-1
h=size(mergezb,1);

% Looks at chromosome and nucleotide posn, looks for variant reads greater
% than 1 percent, for any of the samples in the row, flag as variant calls
% (sites of possibly important nt variation), looks for most common
% nucleotide -> refnt, second most common is variant

% Rather than figure out line by line, import reference human genome seq
% and determine nt posn line by line, how many variants, off by how much
% from the reference

%% Section workb
%Now modify to calculate call totals
workb=zeros(h,s*20);
for i=1:h;
    for j=1:s;
        for k=1:6;
            workb(i,(j-1)*20+k)=mergezb(i,(j-1)*11+k);
        end;
        workb(i,(j-1)*20+7)=mergezb(i,(j-1)*11+11);
        workb(i,(j-1)*20+17)=mergezb(i,(j-1)*11+11);
        workb(i,(j-1)*20+11)=mergezb(i,(j-1)*11+1);
        workb(i,(j-1)*20+12)=mergezb(i,(j-1)*11+2);
        for k=1:4;
            workb(i,(j-1)*20+12+k)=mergezb(i,(j-1)*11+6+k);
        end;
    end;
end;

%spacer
for j=1:(s*2);
    for i=1:h;
        workb(i,(j-1)*10+8)=workb(i,(j-1)*10+3)+workb(i,(j-1)*10+4)+workb(i,(j-1)*10+5)+workb(i,(j-1)*10+6)+workb(i,(j-1)*10+7);
        if workb(i,(j-1)*10+8)<5; % require at least 5 reads at a position to score
            workb(i,(j-1)*10+9)=2;
        else workb(i,(j-1)*10+9)=max(workb(i,((j-1)*10+3):((j-1)*10+7)))/workb(i,(j-1)*10+8);
        end;
    end;
end;

%columns in workb are now: Af denotes read count of As in forward direction, etc.
% 1   2  3  4  5  6   7       8               9     10  11        ---                                 19    20
%chr nt Af Gf Cf Tf indels total-f-read# max-f-read% 0  chr nt Ar Gr Cr Tr indels total-r-read# max-r-read% 0
%                                        or 2 if total-f-read#<10                    similar

%% Section refSeq
% Create a new table for the reference human genome sequence gathered
% from the UCSC genome browser (this genome is hg-19)
% initindex.txt stores the initial index, and this file was created in
% the python analyzer.py script

filetext = fileread('genome.txt');

result = regexprep(filetext, '\s+','');

ntVals = nt2int(result);

% A = 1
% C = 2
% G = 3
% T = 4
ntVals = double(ntVals);

lenRef = length(ntVals);

% Change it to:
% A = 1
% G = 2
% C = 3
% T = 4
% because later, we have that the columns for the counts of A, G, C, T
% are in that order
for i=1:lenRef;
    if(ntVals(i) == 2);
        ntVals(i) = 3;
    elseif(ntVals(i) == 3);
        ntVals(i) = 2;
    end;
end;

refSeq = zeros(lenRef, 2);

for i=1:lenRef;
    refSeq(i, 1) = initIndex + i - 1;
    refSeq(i, 2) = ntVals(i);
end;

if(exon_coord_flush(1, 2) > exon_coord_flush(2, 2));
    caseTSC1 = 1;
    caseTSC2 = 0;
else
    caseTSC1 = 0;
    caseTSC2 = 1;
end;

%% Section noiseRefSeq
% Identify the poly-A sequences and flag them in the reference (run of
% consecutive A's or T's)

noiseRefSeq = [refSeq zeros(lenRef, 1)];

locListA = find(conv(double(ntVals == 1), ones(10, 1), 'valid') == 10);

locListA = transpose(locListA([true diff(locListA) >= 10]));

locListT = find(conv(double(ntVals == 4), ones(10, 1), 'valid') == 10);

locListT = transpose(locListT([true diff(locListT) >= 10]));

locList = sort([locListA; locListT]);

inSeqFlag = 0;

for i=1:lenRef;
    % Get that the letter being repeated is either 1 or 4 (representing A
    % or T respectively)
    if(ismember(i, locList) == 1);
        repeatedVal = refSeq(i, 2);
        inSeqFlag = 1;
    end;

    if(i ~= lenRef);
        if(refSeq(i, 2) == refSeq(i+1, 2) && inSeqFlag == 1);
            noiseRefSeq(i, 3) = 1;
        elseif(inSeqFlag == 1);
            noiseRefSeq(i, 3) = 1;
            inSeqFlag = 0;
        else
            inSeqFlag = 0;
        end;
    end;
end;

%% Section copyWorkb
% Align the reference human genome (hg-19 on the UCSC database) with the nt
% positions considered in the table. The last two columns of the copyWorkb
% table denote the nt and value corresponding to the letter (A, G, C, T)
workbRows = size(workb, 1);

refSeqNew = zeros(workbRows, 2);

copyWorkb = workb;

copyWorkb = [copyWorkb refSeqNew];

[row, col] = find(refSeq == copyWorkb(1, 2));

row = row(1);
col = col(1);

%{
for i=1:workbRows;
    copyWorkb(i, end) = refSeq(row, col+1);
    copyWorkb(i, end-1) = refSeq(row, col);
    row = row + 1;

end;
%}

% If the gene is non-contiguous, then use this (which is slower):

for i=1:workbRows;
    [row, col] = find(refSeq == copyWorkb(i, 2));
    if ~(isempty(row) && isempty(col));
        copyWorkb(i, end) = refSeq(row, col+1);
        copyWorkb(i, end-1) = refSeq(row, col);
    else
    end;
end;


%% Section workNewRef
% Use the reference human genome sequence in order to determine the number
% of variants and by how much this is off from the reference (find the
% percent allele frequency)

% First, set up the new worktable based off of copyWorkb
workNewRef = copyWorkb;

for i=1:s;
    % Six new fields for: A variant allele freq,
    %                     G variant allele freq,
    %                     C variant allele freq,
    %                     T variant allele freq,
    %                     indel variant allele freq,
    %                     max variant allele freq
    newCols = zeros(workbRows, 6);
    workNewRef = [workNewRef(:, 1 : 20 + 20*(i-1) + 6*(i-1)) newCols workNewRef(:, 21 + 20*(i-1) + 6*(i-1) : end)];
end;

workNewRef = [workNewRef zeros(workbRows, s)];

% Next, populate the new fields with the variant allele frequency. This
% table also contains a column for the most common variant nucleotide
% letters in each sample. These columns are at the end of the table

workNewRefWithVals = workNewRef;

for i=1:workbRows;
    % Get the value in the reference genome
    refVal = workNewRef(i, end - s);
    for n=1:s;
        % Get the index of the column corresponding to the ref in the forward dir
        indexFwd = refVal + 2 + (26 * (n - 1));

        % Get the index of the column corresponding to the ref in the reverse dir
        indexRev = indexFwd + 10;

        % Store these indices in an array
        dirVals = [indexFwd indexRev];

        % Get the total number of reads (sum both forward and reverse)
        totReads = workNewRef(i, 8 + 26 * (n - 1)) + workNewRef(i, 18 + 26 * (n - 1));

        if(totReads == 0);
            totReads = 2;
            workNewRef(i, 8 + 26 * (n - 1)) = 1;
            workNewRef(i, 18 + 26 * (n - 1)) = 1;
        end;
        
        % Get the column value of the column corresponding to the reference
        % allele's variant frequency
        refColIndex = 20 + refVal + (26 * (n - 1));

        % Set this to be -1 in order to denote that this was the reference
        % allele
        workNewRefWithVals(i, refColIndex) = -1;

        % Initialize the total number of variants
        varList = zeros(1, 5);

        for j=1:5;
            numVarReads = 0;

            % Get the column value of where we are actually storing the
            % variant frequency
            actualUpdateCol = 20 + j + (26 * (n - 1));

            % Only consider those columns that do not represent the
            % reference genome allele
            if(actualUpdateCol ~= refColIndex);
                % Get the total number of variants for that specific letter
                % (A, G, C, T, or indel)
                
                % Only consider value if there actually was a read at all
                % in each direction
                if(workNewRef(i, 2+j+26*(n-1)) > 0 && workNewRef(i, 12+j+26*(n-1)) > 0);
                    if(j ~= 5);
                        numVarReads = workNewRef(i, j + 2 + (26 * (n - 1))) + workNewRef(i, j + 12 + (26 * (n - 1)));
                    else
                        numVarReads = workNewRef(i, j + 2 + (26 * (n - 1)));
                    end;
                end;

                % Update the place where we are storing the variant allele
                % frequency corresponding to that nt
                workNewRefWithVals(i, actualUpdateCol) = numVarReads / totReads;

                % Update the running sum of the variant frequencies for
                % the current nt
                varList(1, j) = numVarReads / totReads;
            end;
        end;

        % Store the maximum variant allele frequency in the appropriate
        % location
        workNewRefWithVals(i, 26 * n) = max(varList(1, 1:4));

        letterVars = varList(1, 1:4);

        maxVarFreqLetter = find(letterVars==max(max(letterVars)));

        numString = '';

        for ind=1:size(maxVarFreqLetter, 2);
            numString = strcat(numString, int2str(maxVarFreqLetter(1, ind)));
        end;

        maxLetter = str2double(numString);

        if(workNewRefWithVals(i, 26 * n) == 0);
            workNewRefWithVals(i, 26 * s + 2 + n) = 0;
        else
            workNewRefWithVals(i, 26 * s + 2 + n) = maxLetter;
        end;
    end;
end;

%% Section exonLookup
% Create lookup table that contains the data for the exon number and
% distance from nearest exon for all possible nt values from the reference

exonLookup = [workNewRef(:, 1:2) zeros(size(workNewRef, 1), 2)];

for i=1:size(exonLookup, 1);

    numExons = size(T1T2exonsflush, 1);

    foundFlag = 0;
    rowNum = 0;

    index = 0;

    while(index < numExons && foundFlag == 0);
        index = index + 1;

        % Make sure the chr numbers match up
        if(exonLookup(i, 1) == T1T2exonsflush(index, 1));
            if(exonLookup(i, 2) >= T1T2exonsflush(index, 3) && exonLookup(i, 2) <= T1T2exonsflush(index, 4));
                rowNum = index;
                foundFlag = 1;
            end;
        end;
    end;

    % Case when the nucleotide position is in an exonic region
    if(rowNum ~= 0);
        exonLookup(i, 3) = T1T2exonsflush(rowNum, 2);

    % Case when the nucleotide position is not in an exonic region
    else
        % THIS IS FOR SPECIFICALLY THE TSC1 GENE, NEED TO CHANGE THIS FOR
        % TSC2 GENE CASE *************************************************
        distanceVals = zeros(numExons, 2);

        % Put all the distances in a table
        for j=1:numExons;
            distanceVals(j, 1) = abs(T1T2exonsflush(j, 3) - exonLookup(i, 2));
            distanceVals(j, 2) = abs(T1T2exonsflush(j, 4) - exonLookup(i, 2));
        end;

        % Find the absolute value of the distance from the nearest exon
        minElement = min(min(distanceVals));
        [row, col] = find(distanceVals == minElement);

        % Get the actual distance from the nearest exon
        distVal = T1T2exonsflush(row(1), col(1) + 2) - exonLookup(i, 2);

        % Closest exon
        exonLookup(i, 3) = T1T2exonsflush(row(1), 2);

        exonLookup(i, 4) = distVal;
    end;

end;

%% Section addExonData

workNewRefWithVals = [workNewRefWithVals exonLookup(:, 3:4)];

%% Section filternewref
% Filter the data in table workNewRefWithVals based on the minimum variant
% allele frequency and the minimum read count, and both of these values are
% inputted by the user in the python portion of the pipeline

filternewref = zeros(size(workNewRefWithVals, 1), size(workNewRefWithVals, 2));

filterrow = 1;

% Only filter out rows for which at least 1 sample meets the allele
% frequency cutoff and the minimum read count cutoff and have at least 1
% read in both directions

for i=1:workbRows;
    % Flag for whether a row should be ignored by our filtering
    ignoreRow = 0;

    % Number of reads in total for a sample
    totalCounts = 0;
    numunsatisfiedsamples = 0;
    for n=1:s;
        % Get the total number of counts
        forwardcount = workNewRefWithVals(i, 8 + 26 * (n - 1));
        reversecount = workNewRefWithVals(i, 18 + 26 * (n-1));
        totalCounts = forwardcount + reversecount;
        %indelFreq = workNewRefWithVals(i, 25 + 26 * (n - 1));

        % Count the number of samples which do not meet the threshhold
        % readcount and minimum allele frequency or do not have at least 1
        % read in both directions
        if((workNewRefWithVals(i, 26 * n ) < minfreq ) || totalCounts < minreadcount || forwardcount == 0 || reversecount == 0);
            numunsatisfiedsamples = numunsatisfiedsamples + 1;
        else
        end;
    end;

    % If not every sample did not meet the threshhold, then we consider
    % that row
    if(numunsatisfiedsamples ~= s);
        filternewref(filterrow, :) = workNewRefWithVals(i, :);
        filterrow = filterrow + 1;
    else
    end;
end;
%% Section filterwithoutzeros
% Remove rows with zeros from the end of the table
filterwithoutzeros = filternewref(any(filternewref, 2), :);


%% Section restorevarfreqs
numFilterRows = size(filterwithoutzeros, 1);

filterinterm = filterwithoutzeros;

for i=1:numFilterRows;  
    % Get the value in the reference genome
    refVal = filterwithoutzeros(i, end - s - 2);
    
    for n=1:s;
               
        % Get the total number of reads (sum both forward and reverse)
        totReads = filterwithoutzeros(i, 8 + 26 * (n - 1)) + filterwithoutzeros(i, 18 + 26 * (n - 1));

        if(totReads == 0);
            totReads = 1;
        end;
        
        % Get the column value of the column corresponding to the reference
        % allele's variant frequency
        refColIndex = 20 + refVal + (26 * (n - 1));
              
        % Initialize the total number of variants
        varList = zeros(1, 5);
        
        for j=1:5;
            numVarReads = 0;

            % Get the column value of where we are actually storing the
            % variant frequency
            actualUpdateCol = 20 + j + (26 * (n - 1));

            % Only consider those columns that do not represent the
            % reference genome allele
            if(actualUpdateCol ~= refColIndex);
                % Get the total number of variants for that specific letter
                % (A, G, C, T, or indel)
                
                if(j ~= 5);
                    numVarReads = filterwithoutzeros(i, j + 2 + (26 * (n - 1))) + filterwithoutzeros(i, j + 12 + (26 * (n - 1)));
                else
                    numVarReads = filterwithoutzeros(i, j + 2 + (26 * (n - 1)));
                end;
                
                filterinterm(i, actualUpdateCol) = numVarReads / totReads;

                % Update the running sum of the variant frequencies for
                % the current nt
                varList(1, j) = numVarReads / totReads;
            end;
        end;

        % Store the maximum variant allele frequency in the appropriate
        % location
        filterinterm(i, 26 * n) = max(varList(1, 1:4));

          
    end;
end;

%% Section import Genome AD refnt and variant info
genomeADalts = dlmread('reformatgenomeAD.txt');
numLetters = size(genomeADalts, 1);

genomeADposns = dlmread('genomeADposns.txt');
numPosns = size(genomeADposns, 1);
% Use the positions in the above matrix variable as a lookup table in order
% to determine where the sites of interest are as we go through this table
% later

mergedGenomeData = [genomeADposns genomeADalts];


%% Section formatfilter
% The columns in this output data are:
% 1) chromosome number
% 2) nt positions
% 3) 0 (will later be used to represent chr:nt in Excel
% 4) ref nt value (A G C T = 1 2 3 4)
% 5) AF of sample with highest AF
% 6) Variant nt value for sample with highest variant AF
% 7) Number of samples at or greater than the threshold variant AF
% 8) SNP-AF
% 9) Variant AF for Sample 1
% 10) "          for Sample 2
% ... and so on for each of the samples

% First column is the chromosome numbers
formatfilter = zeros(numFilterRows, 8 + s);

% Second column is the nucleotide positions
formatfilter(:, 1:2) = filterwithoutzeros(:, 1:2);

% Fourth column is the reference nucleotice values, where
% A = 1, G = 2, C = 3, T = 4
formatfilter(:, 4) = filterwithoutzeros(:, (26*s) + 2); % Convert these to letters in excel

genomeADfreqs = dlmread('genomeADfreqs.txt');

% Get the allele frequency for the sample with the highest allele frequency
for i=1:numFilterRows;
    maxFreq = filterwithoutzeros(i, 26);
    sampleWithHighest = 1;

    for n=2:s;
        if(filterwithoutzeros(i, 26*n) > maxFreq);
            maxFreq = filterwithoutzeros(i, 26*n);
            sampleWithHighest = n;
        end;
    end;

    % Set the 5th column to be the frequency of the sample
    % with the highest variant allele frequency
    formatfilter(i, 5) = maxFreq;

    % Set the 6th column to be the variant nt value of the sample with the
    % highest allele frequency
    formatfilter(i, 6) = filterwithoutzeros(i, 26*s + 2 + sampleWithHighest);

    numSamplesGeq = 0;

    for n=1:s;
        if(filterinterm(i, 26*n) > 0 && filterinterm(i, 8 + 26 * (n - 1)) + filterinterm(i, 18 + 26 * (n - 1)) >= minreadcount);
            numSamplesGeq = numSamplesGeq + 1;
        else
        end;
    end;

    % Set the 7th column to be the number of samples at or higher than the
    % threshold variant allele frequency
    formatfilter(i, 7) = numSamplesGeq;

    % Set the 8th column to be the percent of variant allele in human
    % population

    [rowAD, colAD] = find(mergedGenomeData == formatfilter(i, 2));

    % Compare to see whether the variant matches the GenomAD database
    % variant changes, and if not, then don't display that value
    displayFlag = 0;
    for index=1:size(rowAD, 1);
        if(mergedGenomeData(rowAD(index, 1), 2) == formatfilter(i, 6));
            displayFlag = 1;
        end;
    end;


    if(displayFlag == 1);
        formatfilter(i, 8) = genomeADfreqs(rowAD(end, 1), 1);
    end;

    % Set the 9th, 10th, and so on columns to be the variant allele
    % frequency for each sample
    for m=1:s;
        formatfilter(i, 8+m) = filterinterm(i, 26*m);
    end;

end;

%% Section formatfilterexons

% Now column 9 represents the exon number at that nucleotide position and
% column 10 represents the distance that this nucleotide position is from
% the nearest exon. Columns 11, 12, etc. represent the variant allele
% frequency for each for the samples. Exon number of 0 means that this
% nucleotide is not in an exonic region.
formatfilterexons = [formatfilter(:, 1:8) zeros(numFilterRows, 2) formatfilter(:, 9:end)];

for i=1:numFilterRows;

    formatfilterexons(i, 9) = filterwithoutzeros(i, end-1);

    formatfilterexons(i, 10) = filterwithoutzeros(i, end);

end;

%% Section nt_coding
filetext2 = fileread('nt_coding.txt');

resultLetters = regexprep(filetext2, '\s+','');

codingNtVals = nt2int(resultLetters);

% A = 1
% C = 2
% G = 3
% T = 4
codingNtVals = double(codingNtVals);

lenCode = length(codingNtVals);

% Change it to:
% A = 1
% G = 2
% C = 3
% T = 4
% because later, we have that the columns for the counts of A, G, C, T
% are in that order
for i=1:lenCode;
    if(codingNtVals(i) == 2);
        codingNtVals(i) = 3;
    elseif(codingNtVals(i) == 3);
        codingNtVals(i) = 2;
    end;
end;

nt_coding = transpose(codingNtVals);

TSC1_nt_coding = nt_coding;
TSC2_nt_coding = nt_coding;

%% Section cdnaLookup
% Here we build the first part of the lookup table for the cdna data

% The columns in the cdnaLookup table are as follows:
% 1) cDNA number
% 2) nt position
% 3) ref nt value
% 4) codon with A substituted at the index
% 5) codon with G substituted at the index
% 6) codon with C substituted at the index
% 7) codon with T substituted at the index

% Set to 1 if we are doing this analysis for TSC2
%caseTSC2 = 0;

indexIntT1 = 1;

while(mod(exon_coord_flush(indexIntT1, 2), 1) ~= 0)
    indexIntT1 = indexIntT1 + 1;
end;

indexIntT2 = size(exon_coord_flush, 1);

while(mod(exon_coord_flush(indexIntT2, 2), 1) ~= 0)
    indexIntT2 = indexIntT2 - 1;
end;

if(caseTSC1 == 1);
    numCdnaRows = TSC1_exon_coord_flush(indexIntT1,4) - TSC1_exon_coord_flush(indexIntT1,3) + TSC1_exon_coord_flush(indexIntT1,5);
    cdnaLookup = zeros(numCdnaRows, 12);
    exonData = TSC1_exon_coord_flush;
    cdnaIndex = numCdnaRows;

elseif(caseTSC2 == 1);
    numCdnaRows = TSC2_exon_coord_flush(indexIntT2,4) - TSC2_exon_coord_flush(indexIntT2,3) + TSC2_exon_coord_flush(indexIntT2,5);
    cdnaLookup = zeros(numCdnaRows, 12);
    exonData = TSC2_exon_coord_flush;
    cdnaIndex = 1;
end;

rowIndex = 1;

for i=1:size(exonData, 1);
    if(exonData(i, 5) ~= 0);
        for val=exonData(i, 3):exonData(i, 4);
            cdnaLookup(rowIndex, 1) = cdnaIndex;
            cdnaLookup(rowIndex, 2) = val;

            if(caseTSC1 == 1);
                % If TSC1, use this:
                cdnaLookup(rowIndex, 3) = TSC1_nt_coding(cdnaIndex, 1);
            elseif(caseTSC2 == 1);
                % If TSC2, use this:
                cdnaLookup(rowIndex, 3) = TSC2_nt_coding(cdnaIndex, 1);
            end;

            aaNum = ceil(cdnaIndex / 3);
            cdnaLookup(rowIndex, 4) = aaNum;

            rowIndex = rowIndex + 1;

            if(caseTSC1 == 1);
                % For TSC1 use this since it is in anti-sense orientation:
                cdnaIndex = cdnaIndex - 1;
            elseif(caseTSC2 == 1);
                % For TSC2 use this since it is in sense orientation:
                cdnaIndex = cdnaIndex + 1;
            end;

        end;
    end;
end;

% FOR TSC2, do this:
% for i=1:size(cdnaLookup, 1);

loopStart = 0;
loopStep = 0;
loopEnd = 0;
checkInd = 0;

if(caseTSC1 == 1);
    loopStart = size(cdnaLookup,1);
    loopStep = -1;
    loopEnd = 1;
    checkInd = 0;
elseif(caseTSC2 == 1);
    loopStart = 1;
    loopStep = 1;
    loopEnd = size(cdnaLookup,1);
    checkInd = 1;
end;

for i=loopStart:loopStep:loopEnd;
    codonInd = mod(i, 3);

    % For TSC2, do codonInd == 1, for TSC1 do codonInd == 0
    if(codonInd == checkInd);

        if(caseTSC1 == 1);
            firstNTval = num2str(cdnaLookup(i, 3));
            secondNTval = num2str(cdnaLookup(i - 1, 3));
            thirdNTval = num2str(cdnaLookup(i - 2, 3));

        elseif(caseTSC2 == 1);
            firstNTval = num2str(cdnaLookup(i, 3));
            secondNTval = num2str(cdnaLookup(i + 1, 3));
            thirdNTval = num2str(cdnaLookup(i + 2, 3));
        end;

        % For TSC2, do i, i+1, i+2 rather than -
    end;

    codonStr = strcat(firstNTval, secondNTval, thirdNTval);

    for j=1:4;
        if(caseTSC1 == 1);

            if(codonInd == 0);
                codonStr(1) = num2str(j);
            elseif(codonInd == 2);
                codonStr(2) = num2str(j);
            elseif(codonInd == 1);
                codonStr(3) = num2str(j);
            end;

        elseif(caseTSC2 == 1);

            %USE FOR TSC2 CASE:

            if(codonInd == 0);
                codonInd = 3;
            end;

            codonStr(codonInd) = num2str(j);

        end;

        cdnaLookup(i, 4+j) = str2double(codonStr);

        % Get the nucleotide values in the codon

        val1 = str2double(codonStr(1));
        val2 = str2double(codonStr(2));
        val3 = str2double(codonStr(3));

        for row=1:size(nt_aa_transl, 1);
            if(nt_aa_transl(row, 1) == val1 && nt_aa_transl(row, 2) == val2 && nt_aa_transl(row, 3) == val3);
                cdnaLookup(i, 8+j) = nt_aa_transl(row, 4);
                break
            end;
        end;
    end;
end;

%% Section varCodons
% Get the variant codons and put them into a table, which we then populate
% with mutaa numbers

% The columns of the varCodons table are as follows:
% 1) ref nt
% 2) variant nt
% 3) nt position
% 4) codon number
% 5) amino acid index number
% 6) ref amino acid (aa)
% 7) mut amino acid (aa)

varCodons = zeros(size(formatfilterexons, 1), 7);

varCodons(:, 1) = formatfilterexons(:, 4);
varCodons(:, 2) = formatfilterexons(:, 6);
varCodons(:, 3) = formatfilterexons(:, 2);

for i=1:size(varCodons, 1);
    [row, col] = find(cdnaLookup == varCodons(i, 3));

    if ~(isempty(row) && isempty(col));
        varCodons(i, 4) = cdnaLookup(row, 1);

        varCodons(i, 5) = cdnaLookup(row, 4);

        colVal = 8 + varCodons(i, 1);
        varCodons(i, 6) = cdnaLookup(row, colVal);

        if(length(num2str(varCodons(i, 2))) == 1);
            mutCol = 8 + varCodons(i, 2);
            varCodons(i, 7) = cdnaLookup(row, mutCol);
        end;
    end;
end;

%% Section formatfilteraa
% Contains the formatted data that is used in the analysis. Here the
% column values are:
% 1) chromosome number
% 2) nt positions
% 3) 0 (will later be used to represent chr:nt in Excel
% 4) ref nt value (A G C T = 1 2 3 4)
% 5) AF of sample with highest AF
% 6) Variant nt value for sample with highest variant AF
% 7) Number of samples with variant AF at or greater than 0
% 8) SNP Allele Frequency
% 9) exon number
% 10) position (distance from nearest exon)
% 11) codon number
% 12) amino acid index number
% 13) ref aa
% 14) mut aa
% 15) is the nt posn +/- 5 bp away from a polyA sequence?
% 16) Variant AF for Sample 1
% 17) "          for Sample 2
% ... and so on for each of the samples

formatfilteraa = [formatfilterexons(:, 1:10) varCodons(:, 4:7) zeros(size(varCodons, 1), 1) formatfilterexons(:, 11:end)];

for i=1:size(varCodons, 1);
    [row, col] = find(formatfilteraa(i, 2) == noiseRefSeq);

    if(row <= lenRef - 5);
        if(noiseRefSeq(row - 5, 3) == 1 ...
               || noiseRefSeq(row - 4, 3) == 1 ...
               || noiseRefSeq(row - 3, 3) == 1 ...
               || noiseRefSeq(row - 2, 3) == 1 ...
               || noiseRefSeq(row - 1, 3) == 1 ...
               || noiseRefSeq(row, 3) == 1 ...
               || noiseRefSeq(row + 1, 3) == 1 ...
               || noiseRefSeq(row + 2, 3) == 1 ...
               || noiseRefSeq(row + 3, 3) == 1 ...
               || noiseRefSeq(row + 4, 3) == 1 ...
               || noiseRefSeq(row + 5, 3) == 1);

           formatfilteraa(i, 15) = 1;
        end;
    end;
end;

snvdata = [formatfilteraa filterwithoutzeros(:, end-s-1:end-2)];

%% Section exportformatfilter
% Here, we export the data in the formatfilter table to an excel document

excelfile = 'snvdata.csv';
dlmwrite(excelfile, snvdata, 'precision', 9);

%% Section workNewRefIndels

workNewRefIndels = [workNewRefWithVals(:, 1) workNewRefWithVals(:, end-s-3:end-s-2) zeros(workbRows, s * 2) workNewRefWithVals(:,end-1:end)];

for i=1:workbRows;
    filterFlag = 0;
    for n=1:s;
        if(workNewRefWithVals(i, 7 + (26 * (n - 1))) >= 2 && workNewRefWithVals(i, 25 + (26 * (n - 1))) >= minIndelFreq);
            filterFlag = 1;
            break;
        end;
    end;

    if(filterFlag == 1);
        for m=1:s;
            % Here the indel frequency is the maximum of the frequencies of
            % the indels in the two reads (forward and reverse)
            workNewRefIndels(i,3+m) = workNewRefWithVals(i, 7 + (26 * (m - 1)));
            %workNewRefIndels(i,3+m+s) = max((workNewRefWithVals(i, 7 + (26 * (m - 1)))/(workNewRefWithVals(i, 8 + (26 * (m - 1))))), (workNewRefWithVals(i, 17 + (26 * (m - 1)))/(workNewRefWithVals(i, 18 + (26 * (m - 1))))));
            totalReads = ((workNewRefWithVals(i, 8 + (26 * (m - 1)))) + (workNewRefWithVals(i, 18 + (26 * (m - 1)))));
            
            if(totalReads == 0);
                totalReads = 1;
            end;
            
            workNewRefIndels(i, 3+m+s) = (workNewRefWithVals(i, 7 + (26 * (m - 1)))) / totalReads;
        end;
    end;
end;

%% Section filterIndels

filterIndels = zeros(size(workNewRefIndels, 1), size(workNewRefIndels, 2));
index = 1;
for i=1:workbRows;
    keepFlag = 0;

    for n=1:s;
        if(workNewRefIndels(i, 3 + n) >= 2 && workNewRefIndels(i, 3 + n + s) >= minIndelFreq && (((workNewRefWithVals(i, 8 + (26 * (n - 1)))) + (workNewRefWithVals(i, 18 + (26 * (n - 1))))) >= minreadcount));
            keepFlag = 1;
            break;
        end;
    end;

    if(keepFlag == 1);
        filterIndels(index, :) = workNewRefIndels(i, :);
        index = index + 1;

    end;

end;

filterIndelsWithoutZeros = filterIndels(any(filterIndels, 2), :);

%% Section formatfilterindels
% The columns in this formatted output are:
% 1) chromosome
% 2) nt
% 3) 0 (but an excel macro will populate it with chr:nt)
% 4) refnt
% 5) max indel freq
% 6) num samples >= cutoff indel freq
% 7) Exon number (or nearest exon)
% 8) Position from nearest exon
% 9) cDNA number (c#)
% 10) amino acid index number (aa#)
% 11) refnt letter (A G C T = 1 2 3 4)
% 12) Num indels for sample 1
% 13) Num indels for sample 2 ...
% 14) Indel freq for sample 1
% 15) Indel freq for sample 2 ...
% Last two columns: exon number, distance from nearest exon


indelRows = size(filterIndelsWithoutZeros, 1);

formatfilterindels = [filterIndelsWithoutZeros(:, 1:2) zeros(indelRows, 1) filterIndelsWithoutZeros(:, 3) zeros(indelRows, 7) filterIndelsWithoutZeros(:, 4:end)];

for i=1:indelRows;
    formatfilterindels(i, 5) = max(formatfilterindels(i, 12+s:end-2));

    numSamplesAboveCutoff = 0;

    for n=1:s;
        if(formatfilterindels(i, 11+n) >= 2 && formatfilterindels(i, 11+n+s) >= minIndelFreq);
            numSamplesAboveCutoff = numSamplesAboveCutoff + 1;
        end;
    end;

    formatfilterindels(i, 6) = numSamplesAboveCutoff;

    formatfilterindels(i, 7) = filterIndelsWithoutZeros(i, end-1);

    formatfilterindels(i, 8) = filterIndelsWithoutZeros(i, end);

    [row, col] = find(cdnaLookup == formatfilterindels(i, 2));

    if ~(isempty(row) && isempty(col));
        formatfilterindels(i, 9) = cdnaLookup(row, 1);

        formatfilterindels(i, 10) = cdnaLookup(row, 4);

    end;

end;

%% Section indelsPolyA
% Here we display to the user whether the indel took place at or near a
% polyA sequence (+/- 5 bp away from a polyA)

% Now the 12th column of the table represents whether or not tht nt
% position is near a poly-A sequence (1 means that the position is near
% (+/- 5 bp)) a poly-A sequence
indelsPolyA = [formatfilterindels(:, 1:(size(formatfilterindels, 2) - 2 - 2*s)) zeros(indelRows, 1) formatfilterindels(:, (size(formatfilterindels, 2) - 2*s -2 + 1):end-2)];

for i=1:indelRows;
   [row, col] = find(formatfilterindels(i, 2) == refSeq);

   if(row <= lenRef - 5);
       if(noiseRefSeq(row - 5, 3) == 1 ...
               || noiseRefSeq(row - 4, 3) == 1 ...
               || noiseRefSeq(row - 3, 3) == 1 ...
               || noiseRefSeq(row - 2, 3) == 1 ...
               || noiseRefSeq(row - 1, 3) == 1 ...
               || noiseRefSeq(row, 3) == 1 ...
               || noiseRefSeq(row + 1, 3) == 1 ...
               || noiseRefSeq(row + 2, 3) == 1 ...
               || noiseRefSeq(row + 3, 3) == 1 ...
               || noiseRefSeq(row + 4, 3) == 1 ...
               || noiseRefSeq(row + 5, 3) == 1);

           indelsPolyA(i, 12) = 1;
       end;
   end;
end;

%% Section export formatfilterindels
% Here, we export the data in the formatfilterindels table to an csv document

excelfile2 = 'formatfilterindels.csv';
dlmwrite(excelfile2, indelsPolyA, 'precision', 9);

%%%%%%%%%%%%%%%%%%%%%
%{
%% Section binList
% Creates list of fixed-length bins in TSC1/TSC2,
% Columns of binList are as follows:
% 1) bin index
% 2) start nt for the bin
% 3) end nt for the bin
% 4) zeros (but will later represent the average read depth across this
%    bin for sample 1)
% 5) zeros (but will later represent the average read depth across this
%    bin for sample 2)

binSize = 100;

binListRows = ceil((2138713-2097990)/binSize) + ceil((135820020-135766735)/binSize);

binList = zeros(binListRows, 5+s);

if(caseTSC1 == 1);
    % For TSC1, use this:
    binInd = 1;
    index = 135766735;
    endInd = 135820020;

elseif(caseTSC2 == 1);
    % For TSC2, use this:
    binInd = ceil((135820020-135766735)/binSize) + 1;
    index = 2097990;
    endInd = 2138713;
end;


%{
% For TSC1, use this:
while(index < endInd);

    if(caseTSC1 == 1);
        % For TSC1, use this:
        binList(binInd, 1) = 9;

    elseif(caseTSC2 == 1);
        % For TSC2, use this:
        binList(binInd, 1) = 16;
    end;

    binList(binInd, 2) = binInd;
    binList(binInd, 3) = index;

    if(index + binSize - 1 > endInd);
        binList(binInd, 4) = endInd;
    else
        binList(binInd, 4) = index + binSize - 1;
    end;

    index = index + binSize;
    binInd = binInd + 1;

end;
%}

workNewRefWithVals = workNewRefWithVals(any(workNewRefWithVals, 2), :);

% Find the start and end posns of the bins by stepping through the
% worknewrefwithvals table and getting every hundredth element. For the end
% condition, find the end value and then determine if the hundredth element
% is after this end value. If so, then cutoff there for the last bin.

for i=1:binListRows;
    totReads = zeros(1, s);

    % TSC1 case
    if(caseTSC1 == 1);


        row = row(1);

        for j=binList(i,3):binList(i,4);
            [row, col] = find(workNewRefWithVals == binList(j,3));
            
            for m=1:s;
                totReads(1, m) = totReads(1, m) + workNewRefWithVals(row, 8 + 26*(m-1)) + workNewRefWithVals(row, 18 + 26*(m-1));
            end;
        end;

        if(i == ceil((135820020-135766735)/binSize))
            for n=1:s;
                % For TSC2, use a different number, calculate this later
                binList(i, 5+n) = totReads(1, n) / (binList(ceil((135820020-135766735)/binSize), 4) - binList(ceil((135820020-135766735)/binSize), 3) + 1);

            end;
        else
            for n=1:s;
                binList(i, 5+n) = totReads(1, n) / binSize;
            end;
        end;

        for k=1:24;
            if((binList(i, 3) >= T1T2exonsflush(k, 3) && binList(i, 3) <= T1T2exonsflush(k, 4)) || (binList(i, 4) >= T1T2exonsflush(k, 3) && binList(i, 4) <= T1T2exonsflush(k, 4)) || (binList(i, 3) <= T1T2exonsflush(k, 3) && binList(i, 4) >= T1T2exonsflush(k, 4)));
                binList(i, 5) = T1T2exonsflush(k, 2);
                break;
            end;
        end;

    % TSC2 case
    elseif(caseTSC2 == 1);

        for j=binList(i,3):binList(i,4);
            [row, col] = find(workNewRefWithVals == binList(i,3));

            for m=1:s;
                totReads(1, m) = totReads(1, m) + workNewRefWithVals(row, 8 + 26*(m-1)) + workNewRefWithVals(row, 18 + 26*(m-1));
            end;

        end;

        if(i == binListRows);
            for n=1:s;
                % For TSC2, use a different number, calculate this later
                binList(i, 5+n) = totReads(1, n) / (binList(end, 4) - binList(end, 3) + 1);

            end;
        else
            for n=1:s;
                binList(i, 5+n) = totReads(1, n) / binSize;
            end;
        end;

        for k=25:67;
            if((binList(i, 3) >= T1T2exonsflush(k, 3) && binList(i, 3) <= T1T2exonsflush(k, 4)) || (binList(i, 4) >= T1T2exonsflush(k, 3) && binList(i, 4) <= T1T2exonsflush(k, 4)) || (binList(i, 3) <= T1T2exonsflush(k, 3) && binList(i, 4) >= T1T2exonsflush(k, 4)));
                binList(i, 5) = T1T2exonsflush(k, 2);
                break;
            end;
        end;

    end;
end;

%%%%%%%%%%%%%%%%%%%%%

% Section exportbindata
% Here, we export the data in the binList table to an csv document

excelfile3 = 'readcountdata.csv';
dlmwrite(excelfile3, binList, 'precision', 9);



%}

%%%%%%%%%%%%%%%%%%%%%

%% Section workc
% Now convert the read counts in columns 3-7 and 13-17 to fractions of total read# rather than counts in workc
workc=zeros(h,s*20);
for i=1:h;
    for j=1:(s*2);
        for k=1:2
            workc(i,(j-1)*10+k)=workb(i,(j-1)*10+k);
        end;
        for k=8:10
            workc(i,(j-1)*10+k)=workb(i,(j-1)*10+k);
        end;
        if workb(i,(j-1)*10+8)~=0
            for k=3:7
                workc(i,(j-1)*10+k)=workb(i,(j-1)*10+k)/workb(i,(j-1)*10+8);
            end;
        else
        end;
    end;
end;
%% Section T1T2exons
%bring in exon coordinates for calculation of read counts per exon
%!!!!!use new precise exon file!!  T1T2exonsflush.txt
%columns are chr exon 1st-nt last-nt
%Use file import

%% Section T1T2readcount
%determine read counts for each exon
%first trim file down to exon data, new file workca

workca=zeros(size(workc,1),s*20);
k=0;
for i=1:size(workc,1);
    for m=1:size(T1T2exonsflush,1);
        if (workc(i,2)>=T1T2exonsflush(m,3))&&(workc(i,2)<=T1T2exonsflush(m,4));
            k=k+1;
            for h=1:s*20;
                workca(k,h)=workc(i,h);
            end;
        else
        end;
    end;
end;

%isolate read counts in forward and reverse here
temp=zeros(k,s*2);
for i=1:k;
    temp(i,1:2)=workca(i,1:2);
    for h=1:s;
        temp(i,h*2+1)=workca(i,(h-1)*20+8);
        temp(i,h*2+2)=workca(i,(h-1)*20+18);
    end;
end;

%k = # lines in workca
%now sum up read counts at each position for each exon
%modification here which works
T1T2readcount=zeros(size(T1T2exonsflush,1),4+s);
for m=1:size(T1T2exonsflush,1);
    for h=1:4;
        T1T2readcount(m,h)=T1T2exonsflush(m,h);
    end;
    for i=1:k;
        if (workca(i,2)>=T1T2exonsflush(m,3))&&(workca(i,2)<=T1T2exonsflush(m,4));
            for h=1:s;
                T1T2readcount(m,h+4)=T1T2readcount(m,h+4)+temp(i,h*2+1)+temp(i,h*2+2);
            end;
        else
        end;
    end;
end;

%now divide by exon size to get average read count per nt
for m=1:size(T1T2exonsflush,1);
    for h=5:(s+4);
        T1T2readcount(m,h)=T1T2readcount(m,h)/(T1T2exonsflush(m,4)-T1T2exonsflush(m,3));
    end;
end;

%%save /Users/djk/Desktop/T1T2readcount.txt T1T2readcount -ascii -double -tabs

%% Section workd
%use allele frequency cutoff of 0.1%
%use read # cutoff of 2

AFcutoff=0.001;
h=size(workc,1);
workd=zeros(h,s*20);
k=0;
for i=1:size(workc,1);
    mnm=1;
    for j=1:(s*2);
        if (workc(i,(j-1)*10+8)-(workc(i,(j-1)*10+8)*workc(i,(j-1)*10+9)))>1.8
            mnm=min(mnm,workc(i,(j-1)*10+9));
        else
        end;
    end;
    if mnm < 1-AFcutoff;
        k=k+1;
        for j=1:(s*20);
            workd(k,j)=workc(i,j);
        end;
    else
    end;
end;

%k
% k = number of rows with potential variants, is used in subsequent steps,
% here = 26108
% steps above eliminate all nt positions in which there is no sample with > 10 reads in either direction
%Lana: need to check this

%% Section0 base
base=0;
for i=1:k;
    for h=1:(s*2);
        if workd(i,(h-1)*10+8)<10;
            basef=0;
        else
            for j=3:6;
                if workd(i,(h-1)*10+j)==workd(i,(h-1)*10+9);
                    base=j-2;
                else
                end;
                base2freq=0;
            end;

            for j=3:6;
                if j~=base+2;
                    if (workd(i,(h-1)*10+j)>base2freq)&&(workd(i,(h-1)*10+j)*workd(i,(h-1)*10+8)>1.5);
                        %identify the 2nd most freq called nt, and confirm that it has been called more than once
                        base2=j-2;
                        base2freq=workd(i,(h-1)*10+j);
                    else
                    end;
                else
                end;
                if base2freq>0.005;
                    basef=base*10+base2;
                else basef=base;
                end;
                workd(i,(h-1)*10+10)=basef;
            end;
        end;
    end;
end;
% workd(1,1:40)
%columns in workd are: Affrac denotes read fraction of As in forward direction, etc.
% 1   2    3      4      5      6      7            8           9      10  11 - 20
%chr nt Affrac Gffrac Cffrac Tffrac indelfrac total-f-read# max-f-frac base-call  same for reverse reads
%                                        or 2 if total-f-read#<10
% base-call is a number; 1 2 3 4 = A G C T; 2 digits mean first is most common call, second is second most common.
%these are calculated separately for the forward and the reverse reads (col's 10 and 20, respectively, etc.)
%here require 5 reads in a direction to make a base call.

%% Section worke
worke=zeros(k,10*s);
for i=1:k;
    for h=1:s;
        worke(i,(h-1)*10+1)=workd(i,(h-1)*20+8); %total number forward reads
        worke(i,(h-1)*10+2)=workd(i,(h-1)*20+10); %base call (forward)
        if workd(i,(h-1)*20+10)==0; %If there is no forward call
            worke(i,(h-1)*10+3)=0; %Forward max genotype is zero
        else
            if workd(i,(h-1)*20+10)<5; %If not an indel
                worke(i,(h-1)*10+3)=workd(i,(h-1)*20+workd(i,(h-1)*20+10)+2); %forward max allele freq is
            else %If it is an indel
                worke(i,(h-1)*10+3)=workd(i,(h-1)*20+round((workd(i,(h-1)*20+10)-1)/10)+2); %it is (forward base call-1)/10+2 ????
                worke(i,(h-1)*10+4)=workd(i,(h-1)*20+mod(workd(i,(h-1)*20+10),10)+2); %forward 2nd highest af is forward base call +2
            end;
        end; %same for reverse
        worke(i,(h-1)*10+6)=workd(i,(h-1)*20+18);
        worke(i,(h-1)*10+7)=workd(i,(h-1)*20+20);
        if workd(i,(h-1)*20+20)==0;
            worke(i,(h-1)*10+8)=0;
        else
            if workd(i,(h-1)*20+20)<5;
                worke(i,(h-1)*10+8)=workd(i,(h-1)*20+workd(i,(h-1)*20+20)+12);
            else
                worke(i,(h-1)*10+8)=workd(i,(h-1)*20+round((workd(i,(h-1)*20+20)-1)/10)+12);
                worke(i,(h-1)*10+9)=workd(i,(h-1)*20+mod(workd(i,(h-1)*20+20),10)+12);
            end;
        end;
    end;
end;

%spacer
% worke(1,1:20)
%spacer
for i=1:k;
    for h=1:s;
        worke(i,(h-1)*10+10)=min(worke(i,(h-1)*10+2),worke(i,(h-1)*10+7));
        if (worke(i,(h-1)*10+10)~=0)&&(worke(i,5)~=0)&&(worke(i,(h-1)*10+10)~=worke(i,5));
            worke(i,5)=max(worke(i,5),worke(i,(h-1)*10+10));
            worke(i,5)=max(worke(i,5),10);
        else
            worke(i,5)=max(worke(i,5),worke(i,(h-1)*10+10));
        end;
    end;
end;
%spacer
% worke(k-10,1:20)
%columns in worke are:
%      1         2      3          4      5      6 - 9                         10
%total-f-read# f-GT f-max-AF f-minor-AF       same for reverse reads     min of f and r base-calls
%            0 means no call
%f-max-AF = forward maximum allele fraction; f-minor-AF = forward second-highest allele fraction
%column 10 is a key value since it is the minimum of the genotype of the forward and the reverse read,
%filtering out many bogus calls that are seen in one orientation only.
% GT is a number; 1 2 3 4 5 = A G C T del; 2 digits mean first is most common call, second is second most common.
% column 5 is the maximum base-call over all samples analyzed.
%spacer

%% Section workf
chrnt=zeros(k,2);
for i=1:k;
    chr=workd(i,1);
    nt=workd(i,2);
    for h=1:(2*s-1);
        chr=max(chr,workd(i,h*10+1));
        nt=max(nt,workd(i,h*10+2));
    end;
    chrnt(i,1)=chr;
    chrnt(i,2)=nt;
end;
%above done to capture the chr and nt for the rows in worke

%spacer
workf=zeros(k,s*10);
for i=1:k;
    workf(i,1)=chrnt(i,1);
    workf(i,2)=chrnt(i,2);
    workf(i,3)=worke(i,5);
    for j=1:s;
        for m=1:4;
            workf(i,(j-1)*9+5+m)=worke(i,(j-1)*10+m);
        end;
        for m=1:5;
            workf(i,(j-1)*9+5+4+m)=worke(i,(j-1)*10+5+m);
        end;
        if workf(i,(j-1)*9+5+9)>9;
            workf(i,4)=workf(i,4)+1;
        else
        end;
    end;
end;
% workf(1:10,1:5)
%columns in workf are: chr nt maxGT(col5 from worke) #samples-het-GT  - then the 9 columns for each sample,
%missing col5 = total-f-read# f-GT f-max-AF f-minor-AF [same for reverse reads] min of f and r genotypes

%% Section workg
%spacer
workg=zeros(k,s*10);
j=1;
for i=1:k;
    if workf(i,3)>9;
        for m=1:(s*9+5);
            workg(j,m)=workf(i,m);
        end;
        j=j+1;
    else
    end;
end;
j=j-1;
%j
%j is now the size of the file of genotype calls of interest, having eliminated all retained positions not het
%j = 6752 this time
% workg;

%% Section T1T2SNPs
%open a file of common TSC1/TSC2 SNPs from the exome variant server

%now work on indicator col5 -
%value indicates the allele frequency if a known SNP per 1000g data and per
%exac depending on the file
%T1T2SNPSExonIntron is for 1000 g data and exac compiled
for i=1:j;
    for m=1:size(TSC1TSC2GenomAD,1);
        if workg(i,2)==TSC1TSC2GenomAD(m,2);
            workg(i,5)=TSC1TSC2GenomAD(m,3);
        else
        end;
    end;
end;

%spacer
% j;
%%save /Users/djk/Desktop/workg-point-12-14.txt workg -ascii -double -tabs
%workg is the full file of all variant calls that seem real, not considering many types of artifacts
% contains five initial columns, as above, then 9 columns per sample, as above.

%% Section workga
workga=zeros(j,s+5);
for i=1:j;
    for m=1:5;
        workga(i,m)=workg(i,m);
    end;
    for jj=1:s;
        if workg(i,5+(jj-1)*9+2)>9;
            if workg(i,5+(jj-1)*9+2)==workg(i,5+(jj-1)*9+6);
                workga(i,5+jj)=(workg(i,5+(jj-1)*9+1)*workg(i,5+(jj-1)*9+4)+workg(i,5+(jj-1)*9+5)*workg(i,5+(jj-1)*9+8))/(workg(i,5+(jj-1)*9+1)+workg(i,5+(jj-1)*9+5));
            else
                if (round((workg(i,5+(jj-1)*9+2)-1)/10)~=mod(workg(i,5+(jj-1)*9+6),10))||(round((workg(i,5+(jj-1)*9+6)-1)/10)~=mod(workg(i,5+(jj-1)*9+2),10));
                    workga(i,5+jj)=-1;
                else
                    workga(i,5+jj)=(workg(i,5+(jj-1)*9+1)*workg(i,5+(jj-1)*9+3)+workg(i,5+(jj-1)*9+5)*workg(i,5+(jj-1)*9+8))/(workg(i,5+(jj-1)*9+1)+workg(i,5+(jj-1)*9+5));
                end;
            end;
        else
        end;
    end;
end;
%%save /Users/djk/Desktop/workga-point-12-14.txt workga -ascii -double -tabs
%workga is condensed down from workg to have just one col per sample with AF for the second most common allele,
%following the initial 5 fields of chr nt maxGT #-het-GT exon-SNP-AF-info
%when the variant alleles in the two directions are discordant, the AF value is made -1

%% Section workgb
%captures some additional information from workga

%initial 7 fields of chr nt maxGT #-het-GT exon-SNP-AF-info max-var-AF
%#-1-values then one col per sample with AF for the most common allele
workgb=zeros(j,s+7);
for i=1:j;
    for m=1:5;
        workgb(i,m)=workga(i,m);
    end;
    for m=8:(s+7);
        workgb(i,m)=workga(i,m-2);
    end;
    for jj=1:s;
        if workgb(i,(jj+7))==-1;
            workgb(i,7)=workgb(i,7)+1;
        else
        end;
    workgb(i,6)=max(workgb(i,6),workgb(i,(jj+7)));
    end;
end;

%% Section workgc
%captures exon and c.# info, intron relative position, aa# and aa change
%first 13 fields are
%chr nt maxGT #-het-GT SNP-AF-info max-minor-AF
% 1  2    3       4       5             6
%#-1.1-values exon/intron position-rel-exon c.#
%      7           8         9               10
%aa# WT-aa   mtn-aa
%11    12      13
%then one col per sample with AF for the most common allele


workgc=zeros(j,s+13);
for i=1:j;
    for m=1:7;
        workgc(i,m)=workgb(i,m);
    end;
    for m=14:(s+13);
        workgc(i,m)=workgb(i,m-6);
    end;
end;
%!!!!!!!!! Import TSC1_exon_coord_flush.txt and TSC2_exon_coord_flush.txt
%import nt_aa_transl.txt, TSC2_nt_coding.txt, TSC1_nt_coding.txt

% TSC1 case
for i=1:j;
    if workgc(i,1)==9;
        flag=0;
        for m=1:size(TSC1_exon_coord_flush,1);
            if (workgc(i,2)>=TSC1_exon_coord_flush(m,3))&&(workgc(i,2)<=TSC1_exon_coord_flush(m,4));
                flag=1; % Flag since this nt is in an exonic region
                workgc(i,8)=TSC1_exon_coord_flush(m,2);
                if (m>1&&m<23);
                    workgc(i,10)=TSC1_exon_coord_flush(m,5)+TSC1_exon_coord_flush(m,4)-workgc(i,2);
                    workgc(i,11)=(workgc(i,10)+2-rem(workgc(i,10)+2,3))/3;
                    tnt1=TSC1_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3));
                    tnt2=TSC1_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3)+1);
                    tnt3=TSC1_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3)+2);
                    for jj=1:size(nt_aa_transl,1);
                        if (tnt1==nt_aa_transl(jj,1))&&(tnt2==nt_aa_transl(jj,2))&&(tnt3==nt_aa_transl(jj,3));
                            workgc(i,12)=nt_aa_transl(jj,4);
                        else
                        end;
                    end;
                    nttemp=TSC1_nt_coding(workgc(i,10));
                    % this is special for TSC1, since it is in anti-sense
                    % orientation
                    nttemp2=5-nttemp;

                    % Get the var nt value and see whether this value is
                    % the complement basepair (complement of A is T and the
                    % complement of G is C)
                    if nttemp2==rem(workgc(i,3),10);
                        ntmtn=workgc(i,3)/10-(rem(workgc(i,3),10)/10);
                    else
                        ntmtn=rem(workgc(i,3),10);
                    end;

                    % Complement basepair calculation
                    ntmtn2=5-ntmtn;

                    % Check whether the position of the nucleotide is the
                    % start of a codon
                    if rem(workgc(i,10)-1,3)==0;
                        tnt1=ntmtn2;

                    else
                        % If the nucleotide position is in the middle of codon,
                        % then do the following
                        if rem(workgc(i,10)-1,3)==1;
                            tnt2=ntmtn2;

                        % If the nucleotide position is at the end of
                        % codon, then do the following
                        else
                            tnt3=ntmtn2;
                        end;
                    end;

                    % Use the look-up matrix to determine what protein
                    % number the codon is encoding
                    for jj=1:size(nt_aa_transl,1);
                        if (tnt1==nt_aa_transl(jj,1))&&(tnt2==nt_aa_transl(jj,2))&&(tnt3==nt_aa_transl(jj,3));
                            workgc(i,13)=nt_aa_transl(jj,4);
                        else
                        end;
                    end;
                else
                end;
            else
            end;
        end;

        % If the nt was not in an exonic region...
        if flag==0;
            e3=0;
            e5=0;
            for m=1:(size(TSC1_exon_coord_flush,1)-1);
                if (workgc(i,2)>TSC1_exon_coord_flush(m,4))&&(workgc(i,2)<TSC1_exon_coord_flush(m+1,3));
                    e5=TSC1_exon_coord_flush(m,4)-workgc(i,2);
                    e3=TSC1_exon_coord_flush(m+1,3)-workgc(i,2);
                    workgc(i,8)=TSC1_exon_coord_flush(m+1,2);
                else
                end;
            end;
            if e3>0;
                if abs(e5)>e3;
                    workgc(i,9)=e3;
                else
                    workgc(i,9)=e5;
                end;
            else
            end;
        else
        end;
    else
    end;
end;

% TSC2 Case
for i=1:j;
    if workgc(i,1)==16;
        flag=0;
        for m=1:size(TSC2_exon_coord_flush,1);
            if (workgc(i,2)>=TSC2_exon_coord_flush(m,3))&&(workgc(i,2)<=TSC2_exon_coord_flush(m,4));
                flag=1;
                workgc(i,8)=TSC2_exon_coord_flush(m,2);
                if (m>1&&m<43);
                    workgc(i,10)=TSC2_exon_coord_flush(m,5)-TSC2_exon_coord_flush(m,3)+workgc(i,2);
                    workgc(i,11)=(workgc(i,10)+2-rem(workgc(i,10)+2,3))/3;
                    tnt1=TSC2_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3));
                    tnt2=TSC2_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3)+1);
                    tnt3=TSC2_nt_coding(workgc(i,10)-rem(workgc(i,10)-1,3)+2);
                    for jj=1:size(nt_aa_transl,1);
                        if (tnt1==nt_aa_transl(jj,1))&&(tnt2==nt_aa_transl(jj,2))&&(tnt3==nt_aa_transl(jj,3));
                            workgc(i,12)=nt_aa_transl(jj,4);
                        else
                        end;
                    end;
                    nttemp=TSC2_nt_coding(workgc(i,10));
                    if nttemp==rem(workgc(i,3),10);
                        ntmtn=workgc(i,3)/10-(rem(workgc(i,3),10)/10);
                    else
                        ntmtn=rem(workgc(i,3),10);
                    end;
                    if rem(workgc(i,10)-1,3)==0;
                        tnt1=ntmtn;
                    else
                        if rem(workgc(i,10)-1,3)==1;
                            tnt2=ntmtn;
                        else
                            tnt3=ntmtn;
                        end;
                    end;
                    for jj=1:size(nt_aa_transl,1);
                        if (tnt1==nt_aa_transl(jj,1))&&(tnt2==nt_aa_transl(jj,2))&&(tnt3==nt_aa_transl(jj,3));
                            workgc(i,13)=nt_aa_transl(jj,4);
                        else
                        end;
                    end;
                else
                end;
            else
            end;
        end;
        if flag==0;
            e3=0;
            e5=0;
            for m=1:(size(TSC2_exon_coord_flush,1)-1);
                if (workgc(i,2)>TSC2_exon_coord_flush(m,4))&&(workgc(i,2)<TSC2_exon_coord_flush(m+1,3));
                    e3=workgc(i,2)-TSC2_exon_coord_flush(m,4);
                    e5=workgc(i,2)-TSC2_exon_coord_flush(m+1,3);
                    workgc(i,8)=TSC2_exon_coord_flush(m,2);
                else
                end;
            end;
            if e3>0;
                if abs(e5)>e3;
                    workgc(i,9)=e3;
                else
                    workgc(i,9)=e5;
                end;
            else
            end;
        else
        end;
    else
    end;
end;

%%save /Users/djk/Desktop/workgc-point-3-15.txt workgc -ascii -double -tabs

%trim down to +-25 of an exon boundary
%optional
workgd=zeros(j,s+13);
h=0;
for i=1:j;
    if (workgc(i,9)<=25)&&(workgc(i,9)>=-25)
        h=h+1;
        workgd(h,1:s+13)=workgc(i,1:s+13);
    else
    end;
end;
%h

%elminate values < 0.01 or -1
%optional
workgct=zeros(size(workgc,1),s+13);
h=0;
for i=1:size(workgc,1);
    if workgc(i,6)>=0.01;
        h=h+1;
        workgct(h,1:s+13)=workgc(i,1:s+13);
        for j=14:s+13;
        if workgct(h,j)<0.01
            workgct(h,j)=NaN;
        else
        end;
        end;
    else
    end;
end;

%%save /Users/djk/Desktop/workgct-point-3-15.txt workgct -ascii -double -tabs

%%save /Users/djk/Desktop/workgd-point-3-15.txt workgd -ascii -double -tabs

%% Indel finding section; files are del+suffix
% Convert down to minimal size to capture indels, 5 columns per sample
% in delc and deld:
% chr nt #total-reads(f+r) #indels fr-indels
% capture max indel freq among all samples in position s*5+1;
%in deld #>0.001 freq in s*5+2, #>0,<0.001 in s*5+3 posn
% deld is trimmed down to only samples with max-indel-freq > 0.001
minindelfreq=0.005; %this is modulable to desired level of detection
minindelcount=2; %this is modulable to desired level of detection
h=size(workb,1);
delc=zeros(h,s*5+1);
deld=zeros(h,s*5+3);
k=0;
for i=1:h;
    maxdelcount=0;
    for j=1:s;
        delc(i,(j-1)*5+1)=workb(i,1);
        delc(i,(j-1)*5+2)=workb(i,2);
        delc(i,(j-1)*5+3)=workb(i,(j-1)*20+8)+workb(i,(j-1)*20+18);
        delc(i,(j-1)*5+4)=workb(i,(j-1)*20+7);
        maxdelcount=max(maxdelcount,delc(i,(j-1)*5+4));
        delc(i,(j-1)*5+5)=delc(i,(j-1)*5+4)/delc(i,(j-1)*5+3);
        delc(i,s*5+1)=max(delc(i,s*5+1),delc(i,(j-1)*5+5));
    end;
    if (delc(i,s*5+1)>minindelfreq&&maxdelcount>minindelcount);
        k=k+1;
        deld(k,1:s*5+1)=delc(i,1:s*5+1);
    else
    end;
end;
%k
%k captures the number of rows (chr:nt posns) for which some sample has
%indel freq > minindelfreq and indel count > 2
%k = 1239 this time


for i=1:k;
    for j=1:s;
        if deld(i,(j-1)*5+5)>minindelfreq;
            deld(i,s*5+2)=deld(i,s*5+2)+1;
        else
            if deld(i,(j-1)*5+5)>0;
                deld(i,s*5+3)=deld(i,s*5+3)+1;
            else
            end;
        end;
    end;
end;

%delep truncates down to indel-freqs only for each sample, after initial 7
%values of chr nt max-indel-freq '#indel-freq>minindelfreq'
%'#indel-freq>0,<minindelfreq' bl bl

delep=zeros(k,s+7);
for i=1:k;
    delep(i,1)=deld(i,1);
    delep(i,2)=deld(i,2);
    delep(i,3)=deld(i,s*5+1);
    delep(i,4)=deld(i,s*5+2);
    delep(i,5)=deld(i,s*5+3);
    for j=1:s;
        delep(i,7+j)=deld(i,(j-1)*5+5);
    end;
end;


%%
% now capture exon/intron posn information, goes into cols 6 and 7

for i=1:k;
    if delep(i,1)==9;
        flag=0;
        for m=1:size(TSC1_exon_coord_flush,1);
            if (delep(i,2)>=TSC1_exon_coord_flush(m,3))&&(delep(i,2)<=TSC1_exon_coord_flush(m,4));
                flag=1;
                delep(i,6)=TSC1_exon_coord_flush(m,2);
            else
            end;
        end;
        if flag==0;
            e3=0;
            e5=0;
            for m=1:(size(TSC1_exon_coord_flush,1)-1);
                if (delep(i,2)>TSC1_exon_coord_flush(m,4))&&(delep(i,2)<TSC1_exon_coord_flush(m+1,3));
                    e5=TSC1_exon_coord_flush(m,4)-delep(i,2);
                    e3=TSC1_exon_coord_flush(m+1,3)-delep(i,2);
                    delep(i,6)=TSC1_exon_coord_flush(m+1,2);
                else
                end;
            end;
            if e3>0;
                if abs(e5)>e3;
                    delep(i,7)=e3;
                else
                    delep(i,7)=e5;
                end;
            else
            end;
        else
        end;
    else
    end;
end;
for i=1:k;
    if delep(i,1)==16;
        flag=0;
        for m=1:size(TSC2_exon_coord_flush,1);
            if (delep(i,2)>=TSC2_exon_coord_flush(m,3))&&(delep(i,2)<=TSC2_exon_coord_flush(m,4));
                flag=1;
                delep(i,6)=TSC2_exon_coord_flush(m,2);
            else
            end;
        end;
        if flag==0;
            e3=0;
            e5=0;
            for m=1:(size(TSC2_exon_coord_flush,1)-1);
                if (delep(i,2)>TSC2_exon_coord_flush(m,4))&&(delep(i,2)<TSC2_exon_coord_flush(m+1,3));
                    e3=delep(i,2)-TSC2_exon_coord_flush(m,4);
                    e5=delep(i,2)-TSC2_exon_coord_flush(m+1,3);
                    delep(i,6)=TSC2_exon_coord_flush(m,2);
                else
                end;
            end;
            if e3>0;
                if abs(e5)>e3;
                    delep(i,7)=e3;
                else
                    delep(i,7)=e5;
                end;
            else
            end;
        else
        end;
    else
    end;
end;

%%save /Users/djk/Desktop/delep-3-15.txt delep -ascii -double -tabs

%delep has 7+s columns:
%chr nt max-indel-freq '#indel-freq>minindelfreq'
%'#indel-freq>0,<minindelfreq' exon relative-position;  then s indel AF

%trim down to +-25 of an exon boundary, and delete indels seen in > s/2.5
%samples

delept=zeros(k,s+7);
h=0;
for i=1:k;
    if (delep(i,7)<=25)&&(delep(i,7)>=-25)&&(delep(i,4)<s/2.5)
        h=h+1;
        delept(h,1:s+7)=delep(i,1:s+7);
    else
    end;
end;
%h

%%save /Users/djk/Desktop/delept-3-15.txt delept -ascii -double -tabs

% %reorder columns in delepr
% %chr nt #indel-AF>0.001 #indel-AF>0,<0.001 exon posn
% %s#-highest-AF s#-2nd-highest s#-3rd-highest indel-AF{1->s samples}
% %chisq p values comparing #1 v #2, #2 v #3, #3 v #4.
% delepr(1:k,1:2)=delep(1:k,1:2);
% delepr(1:k,3:6)=delep(1:k,4:7);
% delepr(1:k,7)=delep(1:k,s+7+1);
% delepr(1:k,8)=delep(1:k,s+7+3);
% delepr(1:k,9)=delep(1:k,s+7+5);
% delepr(1:k,10)=delep(1:k,3);
% delepr(1:k,11:11+s-1)=delep(1:k,8:8+s-1);
% delepr(1:k,s+10+1)=delep(1:k,s+7+2);
% delepr(1:k,s+10+2)=delep(1:k,s+7+4);
% delepr(1:k,s+10+3)=delep(1:k,s+7+6);

% %%save /Users/djk/Desktop/delepr-3-15.txt delepr -ascii -double -tabs

%below is from previous, more complex version, delete for this run
% %use sortrows to figure out s# with highest indel-freq
% %4 columns of freq are s# #total-reads(f+r) #indels fr-indels
% %at end deldp columns are:
% %first s*5+3 are identical to deld; next (s-1)*2 columns are sorted by
% %indel-AF largest to smallest, with the s# and p value of difference in
% %indel-AF between each sample and the next
% %column s*7+4 is a pointer column which starts at 0, and is set to 1 if any
% %of: #indel-reads for each of 1st and 2nd sorted samples is < 3; # samples
% %with indel-AF>0.001 is > 7; or p value difference between 1st
% %and 2nd sorted samples is > 0.01.
%
% deldp=zeros(k,s*7+4);
% deldp(1:k,1:s*5+3)=deld(1:k,1:s*5+3);
% for i=1:k;
%     freq=zeros(s,6);
%     for j=1:s;
%         freq(j,1)=j;
%         freq(j,2)=deldp(i,(j-1)*5+3);
%         freq(j,3)=deldp(i,(j-1)*5+4);
%         freq(j,4)=deldp(i,(j-1)*5+5);
%     end;
%     freqsort=sortrows(freq,-4);
%     for j=1:s-1;
%       E1=(freqsort(j,3)+freqsort(j+1,3))/(freqsort(j,2)+freqsort(j+1,2));
%       freqsort(j,5)=(abs(freqsort(j,3)-E1*freqsort(j,2))-.5)^2+(abs(freqsort(j+1,3)-E1*freqsort(j+1,2))-.5)^2;
%       freqsort(j,6)=1-chi2cdf(freqsort(j,5),1);
%       deldp(i,s*5+3+(j-1)*2+1)=freqsort(j,1);
%       deldp(i,s*5+3+(j-1)*2+2)=freqsort(j,6);
%     end;
%     if max(freqsort(1,3),freqsort(2,3))<3;
%         deldp(i,s*7+4)=1;
%     else
%     end;
%     if deldp(i,s*5+2)>7;  %modulable
%         deldp(i,s*7+4)=1;
%     else
%     end;
% %     if freqsort(1,6)>0.01;  these commands can drop out important
% %     findings
% %         deldp(i,s*7+4)=1;
% %     else
% %     end;
% end;
%
% %now trim out samples for which pointer column value is 1
% deldp2=zeros(k,s*7+4);
% h=0;
% for i=1:k;
%     if deldp(i,s*7+4)==0;
%         h=h+1;
%         deldp2(h,1:s*7+4)=deldp(i,1:s*7+4);
%     else
%     end;
% end;
% h
% k=h;
%
% %delep truncates down to indel-freqs only for each sample, after initial 7
% %values of chr nt max-indel-freq '#indel-freq>0.001' '#indel-freq>0,<0.001'
% %h = 1675 this time
%

% %specific to current run - compress (max) down related samples
% %then count # samples with non-zero values in column 8
% deltr=zeros(k,8);
% for i=1:k;
%     deltr(i,1)=max(delepr(i,11:16));
%     deltr(i,2)=delepr(i,17);
%     deltr(i,3)=max(delepr(i,18:19));
%     deltr(i,4:7)=delepr(i,20:23);
%     for j=1:7;
%         if deltr(i,j)>0;
%             deltr(i,8)=deltr(i,8)+1;
%         else
%         end;
%     end;
% end;
%
% %then add this info to delepr2 as an additional column at end
% %and eliminate rows with value > 2, as likely artifacts, in delepr3
% delepr2(1:k,1:26)=delepr(1:k,1:26);
% delepr2(1:k,27)=deltr(1:k,8);
% m=0;
% for i=1:k;
%     if delepr2(i,27)<3;
%         m=m+1;
%         delepr3(m,1:27)=delepr2(i,1:27);
%     else
%     end;
% end;
% m
% %m=152 in this instance

%To do: integrate new SNPs

% Uncomment the next line when running with the "Run matlab script" button
% in the Python tkinter widget
%quit()
