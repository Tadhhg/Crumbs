function hm_addpath()
%HM_ADDPATH  Put both of this model's code folders on the MATLAB path.
%
%   hm_addpath
%
%   The heating code (GUI_MODEL, HEAT_SOLVE, ...) lives beside this file.  The
%   field code it needs -- SALINE_STOGRYN and MODEL_CONFIG -- lives in a
%   separate folder, because the field model is a project in its own right.
%
%   In the development tree that folder is <root>\matlab, two levels up from
%   here.  In a shared copy it may be anywhere, so this LOOKS FOR IT rather
%   than assuming the layout, and says something useful if it is missing.
%   That is the whole reason this function exists: HEAT_GUI and GUI_MODEL used
%   to hard-code ..\..\matlab, which broke the moment the folder was moved.
%
%   See also HEAT_GUI, GUI_MODEL.

persistent ok
here = fileparts(mfilename('fullpath'));
if ~isempty(ok) && ok && exist('saline_stogryn','file') == 2, return; end

addpath(here);
if exist('saline_stogryn','file') == 2, ok = true; return; end

root = fileparts(here);                         % ...\Heatmodel
up   = fileparts(root);                         % ...\<project root>
cand = { fullfile(up,'matlab')                  % the development tree
         fullfile(root,'matlab')
         fullfile(here,'field')                 % packaged as a subfolder
         fullfile(root,'field')
         fullfile(up,'field')
         root
         up };
for k = 1:numel(cand)
    if exist(fullfile(cand{k},'saline_stogryn.m'),'file') == 2
        addpath(cand{k});  ok = true;  return
    end
end

% Not where it usually is -- search the tree before giving up.
for b = {root, up}
    f = dir(fullfile(b{1},'**','saline_stogryn.m'));
    if ~isempty(f), addpath(f(1).folder);  ok = true;  return; end
end

error('hm_addpath:field', ...
    ['Cannot find the field model.  "saline_stogryn.m" should be in a folder ' ...
     'near\n  %s\nIf you copied these files somewhere, copy the WHOLE folder, ' ...
     'not just the .m files you thought you needed.'], here);
end
