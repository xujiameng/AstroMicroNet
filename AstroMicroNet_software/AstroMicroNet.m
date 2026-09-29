function app=AstroMicroNet(pythonExecutable,interface)
% Launch the integrated historical interface. Existing controls are unchanged.
if nargin<1 || isempty(pythonExecutable)
    pythonExecutable=getenv('ASTROMICRONET_PYTHON');
    if isempty(pythonExecutable),pythonExecutable='python';end
end
if nargin<2,interface='v2';end
assert(~contains(pythonExecutable,'"') && ~contains(pythonExecutable,newline),'Invalid Python path.');
[status,message]=system(sprintf('"%s" -c "import numpy, scipy"',pythonExecutable));
assert(status==0,'AstroMicroNet:Dependencies','Install requirements.txt in the selected Python environment. %s',message);
setenv('ASTROMICRONET_PYTHON',pythonExecutable);
addpath(fileparts(mfilename('fullpath')));
switch lower(interface)
    case 'v2', app=Cal_node_detect_App_v2;
    case 'v1', app=Cal_node_detect_App;
    case 'guide',app=Cal_node_detect;
    otherwise,error('Interface must be v2, v1 or guide.');
end
end
