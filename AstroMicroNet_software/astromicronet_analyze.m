function result = astromicronet_analyze(inputData, options)
% Shared MATLAB entry to the numerical backend included in this package.
% inputData is a zero-diagonal R, or T-by-N processed traces if inputKind='X'.
% options: pythonExecutable, rho, seedKey, gateKey, roiIDs, inputKind,
%          definition, maxLagFrames. The homogeneity screen cannot be disabled.
if nargin < 2, options = struct(); end
if ~isfield(options,'pythonExecutable')
    options.pythonExecutable = getenv('ASTROMICRONET_PYTHON');
    if isempty(options.pythonExecutable), options.pythonExecutable = 'python'; end
end
if ~isfield(options,'rho'), options.rho = 0.20; end
if ~isfield(options,'seedKey'), options.seedKey = 'default'; end
if ~isfield(options,'inputKind'), options.inputKind = 'R'; end
base = tempname;
inputPath = [base '_input.mat']; configPath = [base '_config.json'];
outputPath = [base '_output.mat']; summaryPath = [base '_output.json'];
cleanup = onCleanup(@() removeTemporaryFiles({inputPath,configPath,outputPath,summaryPath}));
payload = struct(); payload.(options.inputKind) = double(inputData);
if isfield(options,'roiIDs'), payload.roi_ids = double(options.roiIDs); end
save(inputPath,'-struct','payload','-v7');
configuration = struct('pipeline',struct('rho',options.rho),'seed_key',options.seedKey);
if isfield(options,'gateKey'), configuration.gate_key = options.gateKey; end
if isfield(options,'definition'), configuration.definition = options.definition; end
if isfield(options,'maxLagFrames'), configuration.max_lag_frames = options.maxLagFrames; end
fid = fopen(configPath,'w','n','UTF-8');
assert(fid>=0,'Cannot create temporary config');
fwrite(fid,jsonencode(configuration),'char'); fclose(fid);
entry = fullfile(fileparts(mfilename('fullpath')),'run_astromicronet.py');
paths = {options.pythonExecutable,entry,inputPath,configPath,outputPath};
for i=1:numel(paths)
    assert(~contains(paths{i},'"') && ~contains(paths{i},newline),'Invalid quote/newline in a command path');
end
command = sprintf('"%s" "%s" --input "%s" --config "%s" --output "%s"',paths{:});
[status,message] = system(command);
assert(status==0,'AstroMicroNet:BackendFailure','%s',message);
result = load(outputPath);
result.final = logical(result.final(:));
result.coreIndices = find(result.final);
result.peripheryIndices = find(~result.final);
result.originalCoreIDs = result.roi_ids(result.final);
result.summary = jsondecode(fileread(summaryPath));
end

function removeTemporaryFiles(paths)
for i=1:numel(paths)
    if isfile(paths{i}), delete(paths{i}); end
end
end
