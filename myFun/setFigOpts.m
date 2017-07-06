function setFigOpts(varargin)
%% Description
%{ 
My personal small function to set my favourite plotting settings
%}

%%
% Default options
options{1} = 14;
options{2} = 13;
options{3} = 1.4;

% Get user-supplied options
if nargin > 0
    options = cell(nargin, 1);
    for ii = 1 : nargin
        options{ii} = varargin{ii};
    end
end
    
%%
set(0,'defaultaxesfontsize', options{1});
set(0,'defaulttextfontsize', options{2});
set(0,'DefaultLineLineWidth', options{3});
set(0,'defaulttextfontweight', 'Bold');
set(0,'DefaultTextInterpreter', 'none');

end