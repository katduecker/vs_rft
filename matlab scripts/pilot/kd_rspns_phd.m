
% Function separates responses into 2 according to photodiode signals
% 1: MISC004 picks up 60 Hz RFT, MISC005 67 Hz
% 2: MISC004 67 Hz, MSC005 60 Hz
function [misc004_rspns, misc005_rspns] = kd_rspns_phd(subj)


% reject weird trials in responses (negative rt)
cfg = [];
specswtr = [];
for b = 1:size(subj.rspns,1)
    % store weird trials
    rt = cell2mat(subj.rspns{b}(2,:));
    
    if find(rt < 0)
        specswtr = [specswtr; b, find(rt < 0)];
    end

end

% replace weird trial specs with zero cell
subj.exp.cltrls = subj.exp.trials;
for w = 1:size(specswtr,1)
    subj.exp.cltrls{1}{specswtr(w,1),specswtr(w,2)} = num2cell(zeros(1,7));
end

% separate response into misc004 60 Hz and misc005 60 Hz

% misc004
% color 1 target and 60 Hz
[rowt colt] = find((cell2mat(cellfun(@(x) x{4} == 1, subj.exp.cltrls{1},'UniformOutput',false)) + cell2mat(cellfun(@(x) x{6} == 60, subj.exp.trials{1},'UniformOutput',false))) == 2);
% color 1 distractor and 60 Hz
[rowd cold] = find((cell2mat(cellfun(@(x) x{5} == 1, subj.exp.cltrls{1},'UniformOutput',false)) + cell2mat(cellfun(@(x) x{7} == 60, subj.exp.trials{1},'UniformOutput',false))) == 2);

b = [rowt;rowd];
trl = [colt;cold];
misc004_rspns = {};
for r = 1:length(b)
    misc004_rspns = [misc004_rspns,subj.rspns{b(r)}(:,trl(r))];
end

% misc005
% color 1 target and 60 Hz
[rowt colt] = find((cell2mat(cellfun(@(x) x{4} == 2, subj.exp.cltrls{1},'UniformOutput',false)) + cell2mat(cellfun(@(x) x{6} == 60, subj.exp.trials{1},'UniformOutput',false))) == 2);
% color 1 distractor and 60 Hz
[rowd cold] = find((cell2mat(cellfun(@(x) x{5} == 2, subj.exp.cltrls{1},'UniformOutput',false)) + cell2mat(cellfun(@(x) x{7} == 60, subj.exp.trials{1},'UniformOutput',false))) == 2);

b = [rowt;rowd];
trl = [colt;cold];
misc005_rspns = {};
for r = 1:length(b)
    misc005_rspns = [misc005_rspns,subj.rspns{b(r)}(:,trl(r))];
end