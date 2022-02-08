%% VS + RFT
% PhD project 2

% convert EDF2MAT and store edf object as structure

% [c] Katharina Duecker      

function el = kd_edf2mat2struct(edfconvpath, subjpth, file_in,file_out)

addpath(edfconvpath)
obj = Edf2Mat(fullfile(subjpth,file_in));

% delete converted file
try delete(fullfile(subjpth,'el_edf2mat.mat'))
catch ME
end
    
props = properties(obj);
for p = 1:numel(props)
    try
    el.(props{p})=obj.(props{p});
    catch ME
    end
end

file_out = fullfile(subjpth,file_out);

if ~exist(file_out)
    save(file_out, "el",'-v7.3')
end