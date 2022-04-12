function varargout  = dpm(varargin)
%DPM Dynamic Programming using Matrix operations version 1.1.1
%
%   OPTIONS = DPM() returns the standard options structure.
%
%   DPM(FUN,NX,NU)  saves a random test model with the filename FUN.m
%   and NX number of states and NU number of inputs.
%
%   [[OUT] DYN] = DPM(FUN,PAR,GRD,PRB,OPTIONS) whith FUN is a function 
%   handle with the model function. The DPM-function requires that the 
%   state and input grids are discretized and limited. The function assumes 
%   equally spaced grids and the GRD structure needed is defined as:
%     GRD.X0{.}    = 'initial state' (only used in forward simulation)
%     GRD.XN{.}.hi = 'final state upper constraint'
%     GRD.XN{.}.lo = 'final state lower constraint'
%     GRD.Nx{.}    = 'number of elements in the state grid'
%     GRD.Xn{.}.hi = 'upper boundary of the state grid'
%     GRD.Xn{.}.lo = 'lower boundary of the state grid'
%     GRD.Nu{.}    = 'number of elements in the input grid'
%     GRD.Un{.}.hi = 'upper boundary of the input grid'
%     GRD.Un{.}.lo = 'lower boundary of the input grid'
%   PRB contains the problem parameters
%     PRB.Ts     = time step
%     PRB.N      = number of time steps in problem (defines the problem 
%                  length)
%     [PRB.N0]   = start time index (only used in forward simulation)
%     [PRB.W{i}] = (optional) vectors with length PRB.T containing
%     disturbances. 
%   PAR is a user defined variable that is sent to the model function.
%     Useful if the model need some parameters defined outside the DPM
%     function.
%   OPTIONS is the options structure containing
%     OPTIONS.Waitbar     = 'off';  (on/off to show or hide waitbars)
%     OPTIONS.Verbose     = 'on'    (on/off command window status
%                                   notification)
%     OPTIONS.Warnings    = 'off';  (on/off to show or hide warnings)
%     OPTIONS.SaveMap     = 'on';   (on/off to save cost-to-go map)
%     OPTIONS.MyInf       = 1e4;    (a big number for infeaible states 
%                                   where out.I=1)
%     OPTIONS.Minimize    = 1;      (one/zero if minimizing or maximizing)
%     OPTIONS.InputType   = 'cd'    (string with the same number of 
%                                    characters as number of inputs.
%                                    contains the character 'c' if input is 
%                                    continuous or 'd' if discrete.
%                                    Default is all continuous.)
%     OPTIONS.BoundaryMethod = 'none', 'Line' or 'LevelSet'
%      if BoundaryMethod = 'Line'
%       OPTIONS.FixedGrid = 0;      (zero/one use grid as defined in GRD or
%                                    adjust the according to boundary line)
%       OPTIONS.Iter      = 10;     (maximum number of iterations when 
%                                    inverting model)
%       OPTIONS.Tol       = 1e-8;   (minimum tolerance when inverting
%                                     model)
%      endif
%     OPTIONS.gN{1}                 Cost matrix at the final time 
%                                    (must be of size(OPTIONS.gN{1}) = 
%                                    [GRD.Nx{1} GRD.Nx{2} ... GRD.Nx{.}]).
%     
%
%   OUT = DPM(DYN,FUN,PAR,GRD,PRB,OPTIONS) only calculates the forward
%   simulation using DYN calculated earlier.
%
%   OUT = DPM(N0,DYN,FUN,PAR,GRD,PRB,OPTIONS) only calculates the forward
%   simulation using DYN calculated earlier but from index N0 
%   (in the range 1->prb.N).
%
%   Up to 5e6 grid points (prod([GRD.Nx{1} GRD.Nx{2} ... GRD.Nx{.} 
%          GRD.Nu{1} GRD.Nu{2} ... GRD.Nu{.}])) works well.
%
%
%   EXAMPLES_______________________________________________________________
%     DPM can be used to get standard options
%        options = dpm();
%     
%     DPM can be used to generate a test model
%        dpm('model_test',3,2);
%   
%   FULL DYNAMIC PROGRAMMING EXAMPLES______________________________________
%   To calculate the optimal control for a test model:
%     dpm('model_test',2,2);
%     grd.Nx{1}    = 11; grd.Xn{1}.lo = -50; grd.Xn{1}.hi = 100; 
%     grd.Nx{2}    = 11; grd.Xn{2}.lo = -50; grd.Xn{2}.hi = 100; 
%     grd.Nu{1}    = 11; grd.Un{1}.lo = -5;  grd.Un{1}.hi = 5; 
%     grd.Nu{2}    = 11; grd.Un{2}.lo = -5;  grd.Un{2}.hi = 5; 
%     
%     grd.X0{1}    = 0;
%     grd.X0{2}    = 0;
%     grd.XN{1}.lo = 0; grd.XN{1}.hi = 20;
%     grd.XN{2}.lo = 0; grd.XN{2}.hi = 20;
%
%     prb.Ts     = 1/5;
%     prb.N      = 200*1/prb.Ts + 1;
%     prb.N0     = 1;
%     prb.W{1}   = rand(1,prb.N);
%     options    = dpm();
%     [out dyn]  = dpm('model_test',[],grd,prb,options);
%
%
%   HINTS FOR DEBUGGING____________________________________________________
%   The model function is called twice before the actual DP algorithm is
%   used in order to determine the number of states, inputs, and outputs. 
%   When adding breakpoints to the model function, please be aware of this.
%
%
%   REFERENCES_____________________________________________________________
%   When using DPM please cite:
%   Sundstrom, O. and Guzzella, L., "A Generic Dynamic Programming Matlab
%   Function", In Proceedings of the 18th IEEE International Conference on 
%   Control Applications, pages 1625-1630, Saint Petersburg, Russia, 2009
%
%   LICENSE________________________________________________________________
%   This Source Code Form is subject to the terms of the Mozilla Public 
%   License, v. 2.0. If a copy of the MPL was not distributed with this 
%   file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%   COPYRIGHT______________________________________________________________
%   Copyright 2008- Institute for Dynamic Systems and Control, Department 
%   of Mechanical and Process Engineering, ETH Zurich
%   Author: Olle L. Sundstrom
%   
%   CHANGELOG______________________________________________________________
%   v.1.0.6 (2010) Useable for MATLAB 2007 -- CD
%   v.1.0.6 (15-Dec-2010): Added OPTIONS.VERBOSE 'on'|{'off'}. Prints
%	status of DPM in the command window. Output equivalent to Waitbar but
%	with significantly less computational effort---SE
%   v.1.0.6 (July 2011) removed use of roundn
%   v.1.1.0 (15-Dec-2011) Added Level-Set-Method as an option and new 
%   forward simulation scheme as an option -- PE
%   v.1.1.1 (11-Oc-2012) fixed error handling in line 2394 --PE
%   v.1.1.2 (19-Juni-2013) added licence agreement header -- PE
try
    
varargout = {};
if nargin == 3
    if ~exist([varargin{1} '.m'],'file')
        str{1}  = '';        
        str{end+1}  = ['function [X, C, I, signals] = ' varargin{1} '(inp,par)'];
        str{end+1}  = '\n%% inp.X{i} states';
        str{end+1}  = '%% inp.U{i} inputs';
        str{end+1}  = '%% inp.W{i} disturbances (as defined in dis-struct)';
        str{end+1}  = '%% inp.Ts   time step';
        str{end+1}  = '%% par      struct including user defined parameters';
%         str{end+1}  = '%% lim.X{i} = ["lower" ; "upper"] (limits of each state)';
        str{end+1}  = '\n%% state update (out.X{i} must be set within model function)';
%         str{end+1}  = 'func = inp.Ts.*(par.a.*(inp.X{1}-inp.X{1}.^2/par.b)-inp.U{1});';
        for i=1:varargin{2}
            x = round(rand*(varargin{2}-1))+1;
            u = round(rand*(varargin{3}-1))+1;
            str{end+1}  = ['X{' num2str(i) '} = ' num2str(rand) '.*(inp.X{' num2str(x) '} + inp.U{' num2str(u) '} + inp.W{1})./inp.Ts + inp.X{' num2str(i) '};'];
        end
        str{end+1}  = '\n%% cost (out.C{1} must be set within model function)';
        str{end+1}  = 'C{1} = -inp.Ts.*inp.U{1};';
        str{end+1}  = '\n%% Infeasibility (out.I [zero=feasible/one=infeasible] must be set within model function)';
        str{end+1}  = '\n%% for example if the state is outside the grid or infeasible combinations of inputs and states occur.';
        str{end+1}  = '\n%% The cost of these state-input combinations will be set to options.MyInf by the DPM function.';
        str{end+1}  = 'I = 0;';
%         for i=2:varargin{2}
%             str{end+1}  = ['I = bitor(out.I,bitor(out.X{' num2str(i) '}<lim.X{' num2str(i) '}.lo, out.X{' num2str(i) '}>lim.X{' num2str(i) '}.hi));'];
%         end
        str{end+1}  = '\n%% store signals (store any other signals in the out struct)';
        for i=1:varargin{3}
            str{end+1} = ['signals.U{' num2str(i) '} = inp.U{' num2str(i) '};'];
        end
        fid = fopen([varargin{1} '.m'], 'w');
        for i=1:length(str)
            fprintf(fid, [str{i} '\n']);
        end
        fclose(fid);  
        return
    else
        error('DPM:Internal','Filename exist in the current directory')
    end
elseif nargin == 4
    model   = varargin{1};
    par     = varargin{2};
    grd     = varargin{3};
    dis     = varargin{4};
    inp     = dpm_get_empty_inp(grd,dis,'zero');
    varargout{1} = dpm_get_empty_out(model,inp,par,grd,'zero');
    return;
elseif nargin == 2
    error('DPM:Internal','Grid generation is not supported anymore, grd = dpm(Xset,Uset);');
%     varargout{1} = dpm_compose_grid(varargin{1},varargin{2});
%     warning('DPM:General','Initial state X0, final state constraints XN, and state limits Xn are set to default values.')
%     return
elseif nargin == 0
    options.Waitbar     = 'off';
    options.Verbose     = 'on';
    options.Warnings    = 'on';
    options.SaveMap     = 'on';
    options.MyInf       = 1e4;
    options.Minimize    = 1;
    options.BoundaryMethod = 'none';
    %options.FixedGrid   = 1;
    %options.Iter        = 10;
    %options.Tol         = 1e-8;
    
    varargout{1} = options;
    return;
elseif nargin == 5
    RunForward  = nargout > 1;
    RunBackward = 1;
    model   = varargin{1};
    par     = varargin{2};
    grd     = varargin{3};
    dis     = varargin{4};
    options = varargin{5};
elseif nargin == 6
    RunForward  = 1;
    RunBackward = 0;
    dyn     = varargin{1};
    model   = varargin{2};
    par     = varargin{3};
    grd     = varargin{4};
    dis     = varargin{5};
    options = varargin{6};
elseif nargin == 7
    RunForward  = 1;
    RunBackward = 0;
    dyn     = varargin{2};
    model   = varargin{3};
    par     = varargin{4};
    grd     = varargin{5};
    dis     = varargin{6};
    options = varargin{7};
    t0      = varargin{1};
end

% Check all inputs
if nargin==4 || nargin==5 || nargin==6 || nargin==7
    % If disturbance vectors are not set
    if ~isfield(dis,'W')
        dis.W = {};
    end
    grd = input_check_grd(grd,dis.N);
end

if exist('options','var')
    if isfield(options,'Warnings') && strcmp(options.Warnings,'on')
        warning('on','DPM:Backward')
        warning('on','DPM:Forward')
        warning('on','DPM:General')
    end
    
    if ~isfield(options,'CalcLine')
        
        %Backward Compatibility:
        if isfield(options,'InfCost')
            warning('DPM:Internal','The option ''InfCost'' has been renamed to ''MyInf''. Consider adjusting your code.')
            options.MyInf = options.InfCost;
            options = rmfield(options,'InfCost');
        end
        if isfield(options,'BoundaryLineMethod')
            warning('DPM:Internal','The option ''BoundaryLineMethod'' has been renamed to ''BoundaryMethod''. Consider adjusting your code.')
            options.BoundaryMethod = options.BoundaryLineMethod;
            options = rmfield(options,'BoundaryLineMethod');
        end
        if isfield(options,'HideWaitbar')
            warning('DPM:Internal','The option ''HideWaitbar'' has been renamed to ''Waitbar''. Consider adjusting your code.')
            if options.HideWaitbar
                options.Waitbar = 'off';
            else
                options.Waitbar = 'on';
            end
            options = rmfield(options,'HideWaitbar');
        end
         if isfield(options,'UseLine')
            warning('DPM:Internal','The option ''UseLine'' has been renamed to ''BoundaryMethod''. Consider adjusting your code.')
            if options.UseLine
                options.BoundaryMethod = 'Line';
            else
                options.BoundaryMethod = 'none';
            end
            options = rmfield(options,'UseLine');
        end
        
        
        % Sanity Check of options structure:
        fnames = fieldnames(options);
        onames = {'Waitbar';'Verbose';'Warnings';'SaveMap';'MyInf';'Minimize';'BoundaryMethod';'Iter';'Tol';'FixedGrid';'gN';'InputType';'CalcLine';'UseUmap';'UseLine';'UseLevelSet'};
        for ii = 1:length(fnames)
            ok = 0;
            for jj=1:length(onames)
                ok = ok || strcmp(fnames{ii},onames{jj});
            end
            if ~ok
                warning('DPM:Internal',['Unknown option: ''',fnames{ii},'''!'])
            end
        end
        clear fnames onames ii jj ok
        
        options.CalcLine = 0;
        
        if ~isfield(options,'SaveMap');
            options.SaveMap = 0;
        else
            %Backwards compatibility
            if ~isnumeric(options.SaveMap) 
                switch options.SaveMap
                    case 'on'
                        options.SaveMap = 1;
                    case 'off'
                        options.SaveMap = 0;
                    otherwise
                        warning('DPM:Internal','Unable to interpret the option "SaveMap". Using SaveMap = ''on''!' )
                        options.SaveMap = 1;
                end
            end
        end
        
        %interpret user input regarding boundary method
        if ~isfield(options,'BoundaryMethod')
            options.BoundaryMethod = '';
        end
        
        switch options.BoundaryMethod
            case 'none'
                options.UseLine = 0;
                options.UseLevelSet = 0;
            case 'Line'
                options.UseLine = 1;
                options.UseLevelSet = 0;
                options.UseUmap = 1;
            case 'LevelSet'
                options.UseLine = 0;
                options.UseLevelSet = 1;
                options.UseUmap = 0;
            otherwise
                disp('Boundary Method not specified! Using BoundaryMethod = ''none''.')
                options.BoundaryMethod = 'none';
                options.UseLine = 0;
                options.UseLevelSet = 0;
        end
        
        if ~isfield(options,'UseUmap')
            options.UseUmap = 1;
        end
        
        % Boundary Line Method
        if options.UseLine
            if length(grd.Nx)>1
                warning('DPM:Internal','Boundary-Line method works only with one-dimensional systems. Consider using the Level-Set method.')
            end
            if ~isfield(options,'FixedGrid')
                warning('DPM:Internal','Grid adaptation unspecified. Using options.FixedGrid=1.')
                options.FixedGrid = 1;
            end
            if ~isfield(options,'Tol')
                warning('DPM:Internal','Model inversion requires the specification of a tolerance. Using options.Tol=10-8.')
                options.Tol = 1e-8;
            end
            if ~isfield(options,'Iter')
                warning('DPM:Internal','Model inversion requires the specification of a maximum number of iterations. Using options.Iter=10.')
                options.Iter = 10;
            end
            
        end
        
        % Level Set Method:
        if options.UseLevelSet
            if length(grd.Nx)==1
                disp('For one-dimensional systems, consider using the Boundary-Line method.')
            end
            if ~options.SaveMap
                warning('DPM:Internal','Level-Set method needs cost-to-go for forward simulation. Setting SaveMap = ''on''.')
                options.SaveMap = 1;
            end
            options.UseUmap = 0;
        end
        
    else
        
    end
    
end


% DYNAMIC PROGRAMMING _____________________________________________________
% Returns the optimal cost-to-go dyn.Jo with the
% size of T x length(X1)
% and the optimal input matrix dyn.Uo with the size
% size of T x length(X1)
if RunBackward
    for i=1:length(grd.Xn)
        for j=1:length(grd.Xn{i})
            if grd.Xn{i}.lo > grd.Xn{i}.hi
                if length(grd.Xn{i}) > 1
                    error('DPM:Internal',['Upper state boundary for state ' num2str(i) ' is at instance ' num2str(j) ' smaller than the lower state boundary.'])
                else
                    error('DPM:Internal',['Upper state boundary for state ' num2str(i) ' is smaller than the lower state boundary.'])
                end
            end
        end
    end
    if ~options.CalcLine
        %grdXN = linspace(grd.Xn{1}.lo(end),grd.Xn{1}.hi(end),grd.Nx{1}(end))';
        
        % remove if 2d boundary line
%         if options.FixedGrid && isempty(find(grdXN>=grd.XN{1}.lo & grdXN<=grd.XN{1}.hi,1))
%             error('DPM:Internal','Final state constraints are too thight. Try to change to FixedGrid=0 \n\t or widen the final constraints or increase state resolution.')
%         else
        if isfield(options,'FixedGrid') && ~options.FixedGrid && min(grd.XN{1}.hi-grd.XN{1}.lo)/eps<grd.Nx{1}(end)
            error('DPM:Internal','Final state constraints are too tight to include grd.Nx{1}(end) points.\n\t Widen the final constraints.')
        end    
    end
    
    dyn = dpm_backward(model,par,grd,dis,options);
end
% EVAULATE RESULT _________________________________________________________
% Uses the optimal input matrix dyn.Uo to simulate
% the optimal trajectory.
% Returns res with all output and inputs of the model
if RunForward
    try
        out = dpm_forward(dyn,model,par,grd,dis,options);
    catch
        clear_waitbars();
        err = lasterror;
        fprintf('DPM:Forward simulation error \n \t Make sure the problem is feasible.\n')
        inp  = dpm_get_empty_inp(grd,dis,'zero');
        out  = dpm_get_empty_out(model,inp,par,grd,'nan');        
%         varargout{length(varargout)+1} = err;
    end
    varargout{length(varargout)+1} = out;
end
varargout{length(varargout)+1} = dyn;
warning('off','DPM:Backward')
warning('off','DPM:General')
warning('off','DPM:Forward')

catch
    err = lasterror;
    notify_user_of_error(err);
    for i=1:nargout
        varargout{i} = [];
    end    
end

function dyn       	= dpm_backward(model,par,grd,dis,options)
% DP_BACKWARD   Computes the optimal input matrix and the optimal cost-to-go.
%
%    Examples:
%
%       Computes the optimal input matrix and the optimal cost-to-go
%          dyn = dpm_backward(model,par,grd,dis,options)
%
%   -THIS IS THE FIRST CORE OF THE DYNAMIC PROGRAMMING-
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström

% Calculate boundaries
%xs = 0;

% Initialize optimal input, and cost maps
%in = 0;
x_sze = ones(size(grd.Nx));
for i=1:length(grd.Nx)
    x_sze(i) = grd.Nx{i}(end);
end
u_sze = ones(size(grd.Nu));
for i=1:length(grd.Nu)
    u_sze(i) = grd.Nu{i}(end);
end

% If UseLine ______________________________________________________________
if options.UseLine
    dyn.B.lo  = dpm_boundary_line_lower(model,par,grd,dis,options);
    dyn.B.hi  = dpm_boundary_line_upper(model,par,grd,dis,options);
end        


% generate current state grid
for i=1:length(grd.Nx)
    if options.UseLine && i==1 && ~options.FixedGrid
        if length(grd.Nx)==1
            current_grd.X{i} = linspace(dyn.B.lo.Xo(end),dyn.B.hi.Xo(end),grd.Nx{i}(end))';
        else
            current_grd.X{i} = linspace(min(dyn.B.lo.Xo(:,i)),max(dyn.B.hi.Xo(:,i)),grd.Nx{i}(i))';
        end
    else
        current_grd.X{i} = linspace(grd.Xn{i}.lo(end),grd.Xn{i}.hi(end),grd.Nx{i}(end))';
        if options.CalcLine
            grd.Xn{1}.lo = grd.Xn_true{1}.lo;
            grd.Xn{1}.hi = grd.Xn_true{1}.hi;
        end
    end
end



% generate current input grid
for i=1:length(grd.Nu)
    current_grd.U{i} = linspace(grd.Un{i}.lo(end),grd.Un{i}.hi(end),grd.Nu{i}(end))';
end

inp0  = dpm_get_empty_inp(grd,dis,'zero');
out0  = dpm_get_empty_out(model,inp0,par,grd,'zero');

% Initialize the outputs dyn.Uo and dyn.Jo
for i=1:length(out0.C)
    % if dyn.Jo is specified in options
    if isfield(options,'gN') && length(options.gN)>= i
        if dpm_sizecmp(options.gN{i},zeros(get_size(current_grd)))
            dyn.Jo{i} = options.gN{i};
        else
            error('DPM:Internal',['options.gN{' num2str(i) '} has incorrect dimesions']);
        end
    % if dyn.Jo is NOT specified in options
    else
        dyn.Jo{i} = zeros(get_size(current_grd));
    end
end

% if not using boundary line set cost of states outside feasible region
% to myinf
if ~options.UseLine && ~options.CalcLine && ~options.UseLevelSet
    for i=1:length(grd.Nx)
        eval(['dyn.Jo{1}(' repmat(':,',1,i-1) 'current_grd.X{i} > grd.XN{i}.hi' repmat(',:',1,length(grd.Nx)-i) ') = options.MyInf;'])
        eval(['dyn.Jo{1}(' repmat(':,',1,i-1) 'current_grd.X{i} < grd.XN{i}.lo' repmat(',:',1,length(grd.Nx)-i) ') = options.MyInf;'])
    end
end


if options.UseLevelSet
    %set up level set function  
    dyn.Jo{end+1} = -inf*zeros(get_size(current_grd));
    %initialize by V_N = h(x_N)
    if length(grd.Nx)>1
        code_fin_cst = '[';
        for i=1:length(grd.Xn)
            code_fin_cst = [code_fin_cst, 'x{', num2str(i),'} '];
        end
        code_fin_cst = [code_fin_cst, '] = ndgrid('];
        for i=1:length(grd.Xn)
            code_fin_cst = [code_fin_cst, '(current_grd.X{', num2str(i),'}),'];
        end
        code_fin_cst = code_fin_cst(1:end-1);
        code_fin_cst = [code_fin_cst, ');'];
    else
        code_fin_cst = 'x{1} = current_grd.X{1};';
    end
    eval(code_fin_cst);
    for i=1:length(x)
        dyn.Jo{end} = max(dyn.Jo{end},max(grd.XN{i}.lo-x{i},x{i}-grd.XN{i}.hi)); 
    end
    if length(grd.Nx)>1
        code_x_grd = '[';
        for i=1:length(grd.Xn)
            code_x_grd = [code_x_grd, 'x', num2str(i),' '];
        end
        code_x_grd = [code_x_grd, '] = ndgrid('];
        for i=1:length(grd.Xn)
            code_x_grd = [code_x_grd, '(1:x_sze(', num2str(i),'))'','];
        end
        code_x_grd = code_x_grd(1:end-1);
        code_x_grd = [code_x_grd, ');'];
    else
        code_x_grd = 'x1 = (1:x_sze(1))'';';
    end
end

% initialization if the entire cost-to-go map should be saved
if options.SaveMap
    V_map = cell(length(dyn.Jo),dis.N+1);
    for i=1:length(dyn.Jo)
        for j=1:dis.N+1
            x_sze_n = [];
            for k=1:length(grd.Nx)
                x_sze_n = [x_sze_n grd.Nx{k}(j)];
            end
            if length(x_sze_n)==1
                x_sze_n = [x_sze_n 1];
            end
            V_map{i,j}          = nan(x_sze_n);
        end
        V_map{i,dis.N+1} = dyn.Jo{i};
    end
end

% initialize the optimal input map(s)
dyn.Uo = cell(length(grd.Nu),dis.N);
for i=1:length(grd.Nu)
    for j=1:dis.N
        x_sze_n = [];
        for k=1:length(grd.Nx)
            x_sze_n = [x_sze_n grd.Nx{k}(j)];
        end
        if length(x_sze_n)==1
            x_sze_n = [x_sze_n 1];
        end
        dyn.Uo{i,j} = nan.*ones(x_sze_n);
    end
end


% GENERATE CODE FOR GENERATING GRID
code_generate_grid = '[';
for i=1:length(grd.Nx)
    code_generate_grid = [code_generate_grid 'inp.X{' num2str(i) '} '];
end
for i=1:length(grd.Nu)
    code_generate_grid = [code_generate_grid 'inp.U{' num2str(i) '} '];
end
code_generate_grid = [code_generate_grid '] = ndgrid('];
for i=1:length(grd.Nx)
    code_generate_grid = [code_generate_grid 'current_grd.X{' num2str(i) '},'];
end
for i=1:length(grd.Nu)
    code_generate_grid = [code_generate_grid 'current_grd.U{' num2str(i) '},'];
end
code_generate_grid = code_generate_grid(1:end-1);
code_generate_grid = [code_generate_grid ');'];


% GENERATE CODE FOR COST-TO-GO INTERPOLATION
eval(['xsize = [' dpm_code('length(current_grd.X{#}) ',1:length(grd.Nx)) '];']);
for k=1:length(dyn.Jo)
    code_cost_to_go_interp{k} = ['cost_to_go{' num2str(k) '} = dpm_interpn('];
    for i=fliplr(find(xsize>1))
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'previous_grd.X{' num2str(i) '},'];
    end
    code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'dyn.Jo{' num2str(k) '},'];
    for i=fliplr(find(xsize>1))
        if options.CalcLine
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'inp.X{' num2str(i) '},'];
        else
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'out.X{' num2str(i) '},'];
        end
    end
    code_cost_to_go_interp{k} = code_cost_to_go_interp{k}(1:end-1);
    code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} ');'];
end
% GENERATE CODE FOR CONVERTING IND TO SUB
code_ind2str = '[';
for i=1:length(grd.Nu)
    code_ind2str = [code_ind2str 'uo' num2str(i) ' '];
end
code_ind2str = [code_ind2str '] = ind2sub(u_sze,ui);'];


% DYNAMIC PROGRAMMING
% Display progres bar of Dynamic Programming backward iteration
if ~isdeployed && strcmp(options.Waitbar,'on')
    if ~options.CalcLine
        h = waitbar(1,'DP running backwards. Please wait...');
    else
        h = waitbar(1,'DP calculating boundary line. Please wait...');
    end
    set(h,'name','DPM:Waitbar');
end

if ~isdeployed && strcmp(options.Verbose,'on')
    if ~options.CalcLine
        fprintf('%s','DP running backwards:     %%');
    else
        fprintf('%s','DP calculating boundary line:     %%');
	end
end


if options.UseLine
    if length(grd.Nx)==1
        isfeas = dyn.B.lo.Xo(:,1) > grd.X0{1} | dyn.B.hi.Xo(:,1) < grd.X0{1};
    elseif length(grd.Nx)==2
        isfeas = dpm_interpn(current_grd.X{2},dyn.B.lo.Xo(:,:,1),grd.X0{2}) > grd.X0{1} | dpm_interpn(current_grd.X{2},dyn.B.hi.Xo(:,:,1),grd.X0{2}) < grd.X0{1};
    end
    if sum(reshape(isfeas,1,numel(isfeas)))~=0
        warning('DPM:Backward','Initial value not feasible!')
    end
end

% flag for warnings (only warned one time)
iswarned = 0;

% start of dp-backward loop
% for n=dis.N-1:-1:1
n = dis.N+1;
while n > 1
    n = n-1;
    
    previous_grd = current_grd;
    x_sze = nan(1,length(grd.Nx));
    u_sze = nan(1,length(grd.Nu));
    for i=1:length(grd.Nx)
        if options.UseLine && i==1 && ~options.FixedGrid
            if length(grd.Nx)==1
                current_grd.X{i} = linspace(dyn.B.lo.Xo(n),dyn.B.hi.Xo(n),grd.Nx{i}(n))';
            else
                current_grd.X{i} = linspace(min(min(dyn.B.lo.Xo(1,:,n:n+1))),max(min(dyn.B.hi.Xo(1,:,n:n+1))),grd.Nx{i}(n))';
            end
        elseif ~options.CalcLine || i~=1
            current_grd.X{i} = linspace(grd.Xn{i}.lo(n),grd.Xn{i}.hi(n),grd.Nx{i}(n))';
        end
        x_sze(i) = grd.Nx{i}(n);
    end
    for i=1:length(grd.Nu)
        current_grd.U{i} = linspace(grd.Un{i}.lo(n),grd.Un{i}.hi(n),grd.Nu{i}(n))';
        u_sze(i) = grd.Nu{i}(n);
    end

    % generate input and state grid
    eval(code_generate_grid);

    if options.CalcLine
        % if calculating boundary line initialize the first input state
        eval(['inp.X{1} = repmat(dyn.Jo{1},[ones(1,length(grd.Nx)) ' dpm_code('length(current_grd.U{#}) ',1:length(grd.Nu)) ']);']);
    end

    % Call model function _________________________________________________
    % generate disturbance
    for w = 1:length(dis.W)
        inp.W{w} = dis.W{w}(n);
    end
    inp.Ts   = dis.Ts;

try
    % call model function
    if isfield(options,'Signals') && ~isempty(options.Signals)
        [out.X out.C out.I signals] = feval(model,inp,par);
    else
        [out.X out.C out.I] = feval(model,inp,par);
    end
    
    if options.UseLevelSet
        %Calculate Level set function:
        if n==dis.N %check analytically
            Vt = -inf*ones(size(out.X{1}));
            for i=1:length(grd.Xn)
                Vt = max(Vt,max(grd.XN{i}.lo-out.X{i},out.X{i}-grd.XN{i}.hi));
            end
        else %check by interpolation
            eval(code_cost_to_go_interp{end});
            Vt = cost_to_go{end};
        end
        Vt(out.I==1) = options.MyInf;
        if length(grd.Nu)>1
            Vt = reshape(Vt,[x_sze prod(u_sze)]);
        end
        [dyn.Jo{end} ub] = min(Vt,[],length(x_sze)+1);

    end
    
    % determine the arc-cost
    for i=1:length(grd.Nx)
        out.I = bitor(out.I,out.X{i}>grd.Xn{i}.hi(n+1));
        out.I = bitor(out.I,out.X{i}<grd.Xn{i}.lo(n+1));
    end
    J  = (out.I==0).*out.C{1} + out.I.* options.MyInf;
    % _____________________________________________________________________



    % Calculate cost for entire grid
    
        if options.UseLine
                if length(grd.Nx)==1
                    cost_to_go{1} = dpm_interpf1mb(previous_grd.X{1}, dyn.Jo{1}, out.X{1}, [dyn.B.lo.Xo(n+1) dyn.B.hi.Xo(n+1)],[dyn.B.lo.Jo(n+1) dyn.B.hi.Jo(n+1)],options.MyInf);
                else
                    cost_to_go{1} = dpm_interpf2mb(previous_grd.X{1},previous_grd.X{2}, dyn.Jo{1}, out.X{1},out.X{2}, [dyn.B.lo.Xo(:,:,n+1);dyn.B.hi.Xo(:,:,n+1)],[dyn.B.lo.Jo(:,:,n+1);dyn.B.hi.Jo(:,:,n+1)],options.MyInf);
                end
        else
            if length(dyn.Jo{1})==1
                cost_to_go{1} = dyn.Jo{1};
            else
                eval(code_cost_to_go_interp{1});
            end
        end
        % total cost = arc-cost + cost-to-go!
        Jt = J + cost_to_go{1};

        if options.Minimize
            Jt(Jt>options.MyInf) = options.MyInf;
            if options.CalcLine
                Jt(Jt<grd.Xn{1}.lo(n)) = options.MyInf;
            end
        else
            Jt(Jt<options.MyInf) = options.MyInf;
            if options.CalcLine
                Jt(Jt>grd.Xn{1}.hi(n)) = options.MyInf;
            end
        end

catch
    err = lasterror;
        if exist('out','var')
            if sum(reshape(isnan(out.I),1,numel(out.I))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable I')
            end
            if sum(reshape(isnan(out.C{1}),1,numel(out.C{1}))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable C')
            end
            if sum(reshape(isnan(out.X{1}),1,numel(out.X{1}))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable X')
            end
            err.message = [err.message ' Error in dpm_backward at n=' num2str(n)];
        end
        rethrow(err);
    end
    if length(grd.Nu)>1
        Jt = reshape(Jt,[x_sze prod(u_sze)]);
    end
    % minimize the cost-to-go
    if options.Minimize
        [Q ui] = min(Jt,[],length(x_sze)+1);
    else
        [Q ui] = max(Jt,[],length(x_sze)+1);
    end
    
    if options.UseLevelSet
        %handle infeasible states:
        Q_inf = J + cost_to_go{1};
        Q_inf = reshape(Q_inf,[x_sze prod(u_sze)]);
        eval(code_x_grd);
        switch length(grd.Nx)
            case 1
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,ub));
            case 2
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,ub));
            case 3
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,x3,ub));
            case 4
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,x3,x4,ub));
        end
        Q(min(Vt,[],length(x_sze)+1)>0) = Qinf(min(Vt,[],length(x_sze)+1)>0);
        ui(min(Vt,[],length(x_sze)+1)>0) = ub(min(Vt,[],length(x_sze)+1)>0);
    end

    if ~options.CalcLine && isfield(options,'Signals') && ~isempty(options.Signals)
        for i=1:length(options.Signals)
            try
                eval(['dyn.' options.Signals{i} '{n} = signals.' options.Signals{i} '(sub2ind([x_sze u_sze],' dpm_code('(1:x_sze(#))''',1:length(x_sze)) ',ui));']);
            catch
                warning('DPM:Backward','options.signals element is not found in model output.')
            end
        end
    end
    
   
    if sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q)
        if options.UseLine % if using line then all points can be infeasible as long as eventually one point is inside boundary lines.
            if dyn.B.hi.Jo(1,1,n)==options.MyInf && dyn.B.lo.Jo(1,1,n)==options.MyInf || sum(reshape(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1),1,numel(out.X{1})))>0 && sum(reshape(inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1))>dyn.B.lo.Xo(n) & inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1))<dyn.B.hi.Xo(n),1,numel(inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1)))))>0
                if isfield(options,'DebugMode') && options.DebugMode
                    fprintf('DPM:Model function error \n \t Entering model function at the instance where the error occured \n\t Check if the entire grid generates infeasible solutions.\n')
                    if isa(model, 'function_handle')
                        eval(['dbstop in ' func2str(model) ' at 1']);
                    else
                        eval(['dbstop in ' model ' at 1']);
                    end
                    n = n+1;
                    continue
                end
                warning('DPM:Backward','No feasible solution Q(i,j,..) = Inf   for all i,j,...')
                break;
            end
        else
            if ~options.CalcLine && isfield(options,'DebugMode') && options.DebugMode
                if isa(model, 'function_handle')
                    eval(['dbstop in ' func2str(model) ' at 1']);
                else
                    eval(['dbstop in ' model ' at 1']);
                end
                n = n+1;
                continue
            end
            warning('DPM:Backward','No feasible solution Q(i,j,..) = Inf   for all i,j,...')
            break;
        end            
    end    
    

    % Update optimal cost dyn.Jo with the minimum cost-to-go Q
    if options.UseLine
        if length(grd.Nx)==1
            below  = current_grd.X{1}<dyn.B.lo.Xo(n);
            above  = current_grd.X{1}>dyn.B.hi.Xo(n);
            inside = current_grd.X{1}>=dyn.B.lo.Xo(n) & current_grd.X{1}<=dyn.B.hi.Xo(n);
        else
            below  = current_grd.X{1}<min(dyn.B.lo.Xo(:,:,n));
            above  = current_grd.X{1}>max(dyn.B.hi.Xo(:,:,n));
            inside = current_grd.X{1}>=min(dyn.B.lo.Xo(:,:,n)) & current_grd.X{1}<=max(dyn.B.hi.Xo(:,:,n));
        end
        dyn.Jo{1} = nan(size(Q));

        eval(['dyn.Jo{1}(:' repmat(',:',1,length(grd.Nx)-1) ')     = Q;']);
        % Single State: if all points are infeasible and some points are between
        % boundaries: use interpolation between boundary-data
        if dyn.B.lo.Jo(n)<options.MyInf && dyn.B.hi.Jo(n)<options.MyInf && length(grd.Nx)==1 && sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q) && sum(reshape(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n),1,numel(current_grd.X{1})))>0
            dyn.Jo{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)) = dpm_interpn([dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)],[dyn.B.lo.Jo(n) dyn.B.hi.Jo(n)],  current_grd.X{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)));
        end

        eval(['dyn.Jo{1}(above' repmat(',:',1,length(grd.Nx)-1) ') = options.MyInf;']);
        eval(['dyn.Jo{1}(below' repmat(',:',1,length(grd.Nx)-1) ') = options.MyInf;']);
    else
        dyn.Jo{1}       = Q;
    end

    eval(code_ind2str);

    for i=2:length(out.C)
        Xi{1}=[];
        for j=2:length(grd.Nx)
            Xi{j} = reshape(out.X{j},[prod(x_sze) prod(u_sze)]);
        end
        ind2 = dpm_sub2ind([prod(x_sze) prod(u_sze)],1:prod(x_sze),ui);

        Ci  = reshape(out.C{i},[prod(x_sze) prod(u_sze)]);
        if length(grd.Nx)>1
            eval(['cost_to_go{2} = dpm_interpn(' dpm_code('current_grd.X{#},',2:length(grd.Nx)) 'dyn.Jo{i}' dpm_code(',out.X{#}',2:length(grd.Nx)) ');']);
            dyn.Jo{i} = (reshape(cost_to_go{2}(ind2),x_sze)~=options.MyInf & reshape(out.I(ind2),x_sze)==0).*(reshape(Ci(ind2),x_sze) + reshape(cost_to_go{2}(ind2),x_sze)) + (reshape(cost_to_go{2}(ind2),x_sze)==options.MyInf | reshape(out.I(ind2),x_sze)~=0).*options.MyInf;
        else
            dyn.Jo{i} = (dyn.Jo{i}~=options.MyInf & out.I(ind2)==0).*(Ci(ind2) + dyn.Jo{i}) + (dyn.Jo{i}==options.MyInf | out.I(ind2)~=0).*options.MyInf;
        end

        if options.Minimize
            dyn.Jo{i}(dyn.Jo{i}>options.MyInf) = options.MyInf;
        else
            dyn.Jo{i}(dyn.Jo{i}<options.MyInf) = options.MyInf;
        end  

    end
    % Store the optimal input that minimized the cost Jt
    if options.UseLine
        for i=1:length(grd.Nu)
            if length(grd.Nx) > 1
                eval(['dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = current_grd.U{i}(uo' num2str(i) ');']);

                [ind col] = dpm_sub2indr([grd.Nx{1}(n) grd.Nx{2}(n)],ones(1,grd.Nx{2}(n)),dpm_findl(current_grd.X{1},dyn.B.lo.Xo(:,:,n)),2);               
                [s1 s2] = ind2sub([grd.Nx{1} grd.Nx{2}],ind);
                dyn.Uo{i,n}(dpm_sub2ind(size(dyn.Uo{i,n}),s1,s2)) = dyn.B.lo.Uo{i}(dpm_sub2ind(size(dyn.B.lo.Uo{i}),ones(size(col)),col,n.*ones(size(col))));

                [ind col] = dpm_sub2indr([grd.Nx{1}(n) grd.Nx{2}(n)],dpm_findu(current_grd.X{1},dyn.B.hi.Xo(:,:,n)),grd.Nx{1}(n).*ones(1,grd.Nx{2}(n)),2);               
                [s1 s2] = ind2sub([grd.Nx{1}(n) grd.Nx{2}(n)],ind);
                dyn.Uo{i,n}(dpm_sub2ind(size(dyn.Uo{i,n}),s1,s2)) = dyn.B.hi.Uo{i}(dpm_sub2ind(size(dyn.B.lo.Uo{i}),ones(size(col)),col,n.*ones(size(col))));
            else


                    % if cost-to-go is infeasible and some grid points are
                    % still inside boundaries
                    eval(['dyn.Uo{i,n}(:) = current_grd.U{i}(uo' num2str(i) ');']);
                    if dyn.B.lo.Jo(n)<options.MyInf && dyn.B.hi.Jo(n)<options.MyInf && sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q) && sum(reshape(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n),1,numel(current_grd.X{1})))>0
                        dyn.Uo{i,n}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)) = dpm_interpn([dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)],[dyn.B.lo.Uo{i}(n) dyn.B.hi.Uo{i}(n)],  current_grd.X{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)));
                    end
                    dyn.Uo{i,n}(below) = dyn.B.lo.Uo{i}(n);
                    dyn.Uo{i,n}(above) = dyn.B.hi.Uo{i}(n);
            end
        end
    else
        for i=1:length(grd.Nu)
            eval(['dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = current_grd.U{i}(uo' num2str(i) ');']);
        end
    end
    
    % store cost-to-go if savemap options is set
    if options.SaveMap
        for i=1:length(dyn.Jo)
            eval(['V_map{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = dyn.Jo{i};']);
        end
    end

    % Update progres bar for the dynamic programming backward
    if ~isdeployed && strcmp(options.Waitbar,'on')
        waitbar(n/dis.N,h);
    end
    if ~isdeployed && strcmp(options.Verbose,'on') && mod(n-1,floor(dis.N/100))==0 && round(100*n/dis.N)<100
        fprintf('%s%2d %%',ones(1,4)*8,round(100*n/dis.N));
    end
end % finish time loop 

% Clear dyn if in low mem mode
if ~options.SaveMap
    dyn.Jo = [];
else
    dyn.Jo = V_map;
end

% Close progres bar
if ~isdeployed && strcmp(options.Waitbar,'on')
    waitbar(1,h)
    close(h)
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('%s Done!\n',ones(1,5)*8)
end

function LineLower  = dpm_boundary_line_lower(model,par,grd,dis,options)

optb               = options;
parb               = par;
grdb               = grd;
grdb.Xn_true{1}.lo = grd.Xn{1}.lo;
grdb.Xn_true{1}.hi = grd.Xn{1}.hi;
grdb.Nx{1}         = ones(1,dis.N+1);
grdb.Xn{1}.lo      = grd.X0{1}.*ones(1,dis.N+1);
grdb.Xn{1}.hi      = grd.X0{1}.*ones(1,dis.N+1);
optb.UseLine       = 0;
optb.SaveMap       = 1;
parb.model         = model;
parb.options.Iter  = options.Iter;
parb.options.Tol   = options.Tol;
optb.CalcLine      = 1;

if isfield(options,'gN')
    warning('options.gN can only be used w/ 1 state boundary')
    optb.gN{1} = grd.XN{1}.lo;%.*sub_x_ones;
    optb.gN{2} = dpm_interpn(grd.X{1},options.gN{1},grdb.X{1});%.*sub_x_ones);
else
    optb.gN{1} = grd.XN{1}.lo;%.*sub_x_ones;
    optb.gN{2} = zeros(size(grd.XN{1}.lo));%zeros(size(sub_x_ones));
end

optb.MyInf        = options.MyInf;
optb.Minimize     = 1; % PEL:2012 - this used to be =options.Minimize; 
optb.Warnings     = 'off';
dynb              = dpm(@dpm_model_inv,parb,grdb,dis,optb);
% convert cellarray to vectors
vsize = size(dynb.Jo);
for j=1:vsize(1)
    V_new{j} = nan(1,length(dynb.Jo{j,1}),vsize(2));
    for i=1:vsize(2)
        V_new{j}(1,:,i) = [dynb.Jo{j,i}];
    end
end
usize = size(dynb.Uo);
for j=1:usize(1)
    O_new{j} = nan(1,length(dynb.Uo{j,1}),usize(2));
    for i=1:usize(2)
%         if i<usize(2)
            O_new{j}(1,:,i) = [dynb.Uo{j,i}];
%         end
    end
end
dynb.Jo = V_new;
dynb.Uo = O_new;
% convert infeasible points to minimum state boundary
for i=1:dis.N
    dynb.Jo{1}(1,dynb.Jo{1}(1,:,i)>grd.Xn{1}.hi(i),i) = grd.Xn{1}.lo(i);
    dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)          = grd.Xn{1}.lo(i).*ones(size(dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)));
end
dynb.Jo{2}(isnan(dynb.Jo{2}) | dynb.Jo{2}>=optb.MyInf) = options.MyInf;
for i=1:length(dynb.Uo)
    dynb.Uo{i}(isnan(dynb.Uo{i})) = 1;
end
% insert lower boundary into original problem definition
LineLower.Xo       = dynb.Jo{1};
LineLower.Uo       = dynb.Uo;
LineLower.Jo       = dynb.Jo{2}; % add final cost term

function LineUpper  = dpm_boundary_line_upper(model,par,grd,dis,options)

optb               = options;
parb               = par;
grdb               = grd;
grdb.Xn_true{1}.lo = grd.Xn{1}.lo;
grdb.Xn_true{1}.hi = grd.Xn{1}.hi;
grdb.Nx{1}         = ones(1,dis.N+1);
grdb.Xn{1}.lo      = grd.X0{1}.*ones(1,dis.N+1);
grdb.Xn{1}.hi      = grd.X0{1}.*ones(1,dis.N+1);
optb.UseLine       = 0;
optb.SaveMap       = 1;
parb.model         = model;
parb.options.Iter  = options.Iter;
parb.options.Tol   = options.Tol;
optb.CalcLine      = 1;

if isfield(options,'gN')
    warning('options.gN can only be used w/ 1 state boundary')
    optb.gN{1} = grd.XN{1}.hi;%.*sub_x_ones;
    optb.gN{2} = dpm_interpn(grd.X{1},options.gN{1},grdb.X{1});%.*sub_x_ones);
else
    optb.gN{1} = grd.XN{1}.hi;%.*sub_x_ones;
    optb.gN{2} = zeros(size(grd.XN{1}.hi));%zeros(size(sub_x_ones));
end
optb.MyInf        = -options.MyInf;
optb.Minimize     = 0; % PEL:2012 this used to be ~options.Minimize; 
optb.Warnings     = 'off';
dynb              = dpm(@dpm_model_inv,parb,grdb,dis,optb);
% convert cellarray to vectors
vsize = size(dynb.Jo);
for j=1:vsize(1)
    V_new{j} = nan(1,length(dynb.Jo{j,1}),vsize(2));
    for i=1:vsize(2)
        V_new{j}(1,:,i) = [dynb.Jo{j,i}];
    end
end
usize = size(dynb.Uo);
for j=1:usize(1)
    O_new{j} = nan(1,length(dynb.Uo{j,1}),usize(2));
    for i=1:usize(2)
%         if i<usize(2)
            O_new{j}(1,:,i) = [dynb.Uo{j,i}];
%         end
    end
end
dynb.Jo = V_new;
dynb.Uo = O_new;
% convert infeasible points to minimum state boundary
for i=1:dis.N
    dynb.Jo{1}(1,dynb.Jo{1}(1,:,i)<grd.Xn{1}.lo(i),i) = grd.Xn{1}.hi(i);
    dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)          = grd.Xn{1}.hi(i).*ones(size(dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)));
end
dynb.Jo{2}(isnan(dynb.Jo{2}) | dynb.Jo{2}<=optb.MyInf) = options.MyInf;
for i=1:length(dynb.Uo)
    dynb.Uo{i}(isnan(dynb.Uo{i}))        = 1;
end
% insert upper boundary into original problem definition
LineUpper.Xo       = dynb.Jo{1};
LineUpper.Uo       = dynb.Uo;
LineUpper.Jo       = dynb.Jo{2}; % add final cost term

function out       	= dpm_forward(dyn,model,par,grd,dis,options)
% DPM_FORWARD   Simulates the model using the optimal input
%
%    Examples:
%
%       Simulates the model using the optimal input
%          out = dpm_forward(dyn,model,par,grd,dis,options)
%
%   -THIS IS THE SECOND CORE OF THE DYNAMIC PROGRAMMING-
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström


if ~isfield(options,'InputType') || length(options.InputType)~= length(grd.Nu)
    options.InputType = repmat('c',1,length(grd.Nu));
end


code_grid_states = '';
for i=length(grd.Nx):-1:1
    code_grid_states = [code_grid_states 'current_grd.X{' num2str(i) '},'];
end
code_input_states = '';
for i=length(grd.Nx):-1:1
    code_input_states = [code_input_states ',inp.X{' num2str(i) '}'];
end
code_states_nearest = 'ixm1';
for i=2:length(grd.Nx)
    code_states_nearest = [code_states_nearest ',ixm' num2str(i)];
end

inp  = dpm_get_empty_inp(grd,dis,'zero');
outn = dpm_get_empty_out(model,inp,par,grd,'zero');

% initialize state
for i=1:length(grd.Nx)
    inp.X{i}  = grd.X0{i};
    outn.X{i} = grd.X0{i};
end


% Display progres bar of Dynamic Programming forward iteration
if ~isdeployed && strcmp(options.Waitbar,'on')
    h = waitbar(0,'DP running forwards. Please wait...');
    set(h,'name','DPM:Waitbar');
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('DP running forwards:     0 %%');
end

% backward compability of dis.N0
if ~isfield(dis,'T0')
    dis.N0 = 1;
end

if ~options.UseUmap
    % GENERATE CODE FOR GENERATING GRID
    if length(grd.Nu)>1
        code_generate_grid = '[';
        for i=1:length(grd.Nu)
            code_generate_grid = [code_generate_grid 'inpt.U{' num2str(i) '} '];
        end
        code_generate_grid = [code_generate_grid '] = ndgrid('];
        for i=1:length(grd.Nu)
            code_generate_grid = [code_generate_grid 'current_grd.U{' num2str(i) '},'];
        end
        code_generate_grid = code_generate_grid(1:end-1);
        code_generate_grid = [code_generate_grid ');'];
    else
        code_generate_grid = 'inpt.U{1} = current_grd.U{1};';
    end
    
    % GENERATE CODE FOR COST-TO-GO INTERPOLATION
    Jsze = size(dyn.Jo);
    for k=1:Jsze(1)
        code_cost_to_go_interp{k} = ['cost_to_go{' num2str(k) '} = dpm_interpn('];
        for i=length(grd.Nx):-1:1
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'next_grd.X{' num2str(i) '},'];
        end
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'dyn.Jo{' num2str(k) ',n+1},'];
        for i=length(grd.Nx):-1:1
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'X{' num2str(i) '},'];
        end
        code_cost_to_go_interp{k} = code_cost_to_go_interp{k}(1:end-1);
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} ');'];
    end
end

iswarned = 0;
% Forward sim
for n=dis.N0:dis.N
    for i=1:length(grd.Nx)
        if options.UseLine && length(grd.Nx)==1 && i==1 && ~options.FixedGrid
            current_grd.X{i} = linspace(dyn.B.lo.Xo(n),dyn.B.hi.Xo(n),grd.Nx{i}(n))';
        elseif options.UseLine && ~options.FixedGrid
            current_grd.X{i} = linspace(min(dyn.B.lo.Xo(:,n)),max(dyn.B.hi.Xo(:,n)),grd.Nx{i}(n))';
        else
            current_grd.X{i} = linspace(grd.Xn{i}.lo(n),grd.Xn{i}.hi(n),grd.Nx{i}(n))';
        end
    end
    for i=1:length(grd.Nu)
        current_grd.U{i} = linspace(grd.Un{i}.lo(n),grd.Un{i}.hi(n),grd.Nu{i}(n))';
    end
    
    %make inp struct
    for w = 1:length(dis.W)
        inp.W{w} = dis.W{w}(n);
    end
    inp.Ts   = dis.Ts;
    
    % if Umap is not used, try all possible u-candidates:
    if ~options.UseUmap
        for i=1:length(grd.Nx)
            next_grd.X{i} = linspace(grd.Xn{i}.lo(n+1),grd.Xn{i}.hi(n+1),grd.Nx{i}(n+1))';
        end
        
        % set up grid of all possible input candidates
        % make inp struct
        for w = 1:length(dis.W)
            inpt.W{w} = dis.W{w}(n);
        end
        inpt.Ts   = dis.Ts;
        eval(code_generate_grid)
        for i=1:length(grd.Nx)
            inpt.X{i} = inp.X{i}*ones(size(inpt.U{1}));
        end
        
        % input to system
        [X C I] = feval(model,inpt,par);
        % take care of bounds
        for i=1:length(grd.Nx)
            I = bitor(I,X{i}>grd.Xn{i}.hi(n+1));
            I = bitor(I,X{i}<grd.Xn{i}.lo(n+1));
        end
        % if line was caluclated take care of bounds again
        if options.UseLine
            I = bitor(I,X{i}>dyn.B.hi.Xo(n+1));
            I = bitor(I,X{i}<dyn.B.lo.Xo(n+1));
        end
        
        %arc cost
        J  = (I==0).*C{1} + I.* options.MyInf;
        
        if options.UseLevelSet
            % get information about feasiblity:
            eval(code_cost_to_go_interp{end});
            J(cost_to_go{end}>0) = options.MyInf;
        end
        
        % minimize total cost
        % Calculate cost for entire grid
        if length(dyn.Jo{1})==1
            cost_to_go{1} = dyn.Jo{1};
        else
            %interpolate from cost to go map
            eval(code_cost_to_go_interp{1});
        end
        Jt = J + cost_to_go{1};
        
        %if no valid input signal is found:
        if options.UseLevelSet && ~any(reshape(cost_to_go{end},numel(cost_to_go{end}),1)<=0 & reshape(I,numel(I),1)==0)
            Jt = cost_to_go{2};
            Jt(I~=0) = options.MyInf;
        end
        
        if options.Minimize
            Jt(Jt>options.MyInf) = options.MyInf;
        else
            Jt(Jt<options.MyInf) = options.MyInf;
        end
        
        if length(grd.Nu)>1
            Jt = reshape(Jt,[1,numel(Jt)]);
        end
        
        % minimize the cost-to-go
        if options.Minimize
            [Q ui] = min(Jt);
        else
            [Q ui] = max(Jt);
        end
        
        % use input that minimizes total cost
        for i=1:length(inpt.U)
            if length(grd.Nu)>1
                inpt.U{i} = reshape(inpt.U{i},[1,numel(inpt.U{i})]);
            end
            inp.U{i} = inpt.U{i}(ui(1));
        end
        
    else
        if options.UseLine
            for i=1:length(grd.Nu)
                xi = cell(1,length(grd.Nx));
                for j=length(grd.Nx):-1:1
                    if options.InputType(i) == 'd' % Discrete
                        if j==1
                            xistr = dpm_code('xi{#},',2:length(current_grd.X));
                            eval(['x1vec = [dyn.B.lo.Xo(1,' xistr 'n); current_grd.X{j}(current_grd.X{j}>dyn.B.lo.Xo(1,' xistr 'n) & current_grd.X{j}<dyn.B.hi.Xo(1,' xistr 'n)); dyn.B.hi.Xo(1,' xistr 'n)];']);
                            [temp xi{j}] = min(abs(x1vec-inp.X{j}));
                            Xin{j} = x1vec(xi{j});
                        else
                            xi{j} = round((inp.X{j}-current_grd.X{j}(1))/(current_grd.X{j}(2)-current_grd.X{j}(1))) + 1;
                            xi{j} = max(xi{j},1);
                            xi{j} = min(xi{j},length(current_grd.X{j}));
                            Xin{j} = current_grd.X{j}(xi{j});
                        end
                    else % Continuous
                        Xin{j} = inp.X{j};
                    end
                end
                
                if length(grd.Nx) > 1
                    inp.U{i}   = dpm_interpf2sbh(current_grd.X{1},current_grd.X{2}, dyn.Uo{i,n}(:,:), Xin{1},Xin{2}, [dyn.B.lo.Xo(:,:,n); dyn.B.hi.Xo(:,:,n)]);
                else
                    inp.U{i}   = dpm_interpf1sbh(current_grd.X{1}, dyn.Uo{i,n}, Xin{1}, [dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)]);
                end
            end
        elseif ~isempty(grd.Nu)
            for i=1:length(grd.Nu)
                if options.InputType(i) == 'c' % Continuous
                    eval(['inp.U{i}   = dpm_interpn(' code_grid_states 'dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ')' code_input_states ');']);
                else % Discrete
                    for j=1:length(grd.Nx)
                        eval(['ixm' num2str(j) ' = round((inp.X{j}-current_grd.X{j}(1))/(current_grd.X{j}(2)-current_grd.X{j}(1)) + 1);']);
                    end
                    eval(['inp.U{i}    = dyn.Uo{i,n}(' code_states_nearest ');']);
                end
            end
        end
    end
    
    %call model with optimal input:
    [X C I outn] = feval(model,inp,par);
    outn.X = inp.X;
    outn.C = C;
    outn.I = I;
    if outn.I~=0 && ~iswarned
        warning('DPM:Forward','Infeasible Solution!')
        iswarned = 1;
    end
    inp.X = X;
    

    if n > dis.N0
        out = dpm_mergestruct(out,outn);
    else
        out = outn;
    end

    % Update progres bar for the Dynamic Programming
    if ~isdeployed && strcmp(options.Waitbar,'on')
        waitbar((n-dis.N0)/(dis.N-dis.N0),h);
    end   
    
    if ~isdeployed && strcmp(options.Verbose,'on') && mod(n-1,floor(dis.N/100))==0 && round(100*n/dis.N)<100
        fprintf('%s%2d %%',ones(1,4)*8,round(100*n/dis.N));
    end

end


for i=1:length(outn.X)
    out.X{i} = [out.X{i} inp.X{i}];
end


% Close progres bar
if ~isdeployed && strcmp(options.Waitbar,'on')
    waitbar(1,h)
    close(h)
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('%s Done!\n',ones(1,5)*8);
end

function y         	= dpm_interpf1mb(xx,yy,A,xlim,ylim,myInf)
%MY_INTERPF1M Computes the 1D interpolation for the given set A
%   using the function YY(xx). USES EXTRAPOLATION
%
%   Y = MY_INTERPF1M(XX,YY,A)
%
%   XX    = Axis vector for the value matrix YY
%   YY    = Value matrix
%   A     = Set of XX-values to be interpolated
%
%   Y     = Set of interpolated values (same structure as A)
%
%   Assumes equally spaced XX
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström

% find grid point just inside lower boundary
Iinl = find(xx>xlim(1),1,'first');
% find grid point just inside upper boundary
Iinu = find(xx<xlim(2),1,'last');

% find interpolation points between lower boundary and closest grid point
Ibel = find(A >=  xlim(1) & A < xx(Iinl));
% find interpolation points between upper boundary and closest grid point
Ibeu = find(A <=  xlim(2) & A > xx(Iinu));
% find interpolation points outside boundary
Iout = find(A <   xlim(1) | A > xlim(2));
% find interpolation points inside boundary
Iin  = find(A >= xx(Iinl) & A <=xx(Iinu));
% find interpolation points between lower and upper boundary
Iinx = find(A >=  xlim(1) & A <=xlim(2));

% initialize output
y = zeros(size(A));

% interpolate as usual with interior points
if ~isempty(Iin)
    y(Iin) = dpm_interpn(xx',yy',A(Iin));
end
% set outside points to inf
if ~isempty(Iout)
    y(Iout)= myInf;
end
% if there are grid points between boundaries
if ~isempty(find(xx<xlim(2) & xx>xlim(1),1))
    % interpolate points between lower boundary and closest feasible grid point
    if ~isempty(Ibel)
        y(Ibel)= dpm_interpn([xlim(1) xx(Iinl)],[ylim(1) yy(Iinl)],A(Ibel));
    end
    % interpolate points between upper boundary and closest feasible grid point
    if ~isempty(Ibeu)
        y(Ibeu)= dpm_interpn([xx(Iinu) xlim(2)],[yy(Iinu) ylim(2)],A(Ibeu));
    end
else
    % if there are no grid points between boundaries
    y(Iinx)= dpm_interpn(xlim,ylim,A(Iinx));
end
function y         	= dpm_interpf2mb(xx1,xx2,YY,A1,A2,xlim,ylim,myInf)
%MY_INTERPF1M Computes the 1D interpolation for the given set A
%   using the function YY(xx). USES EXTRAPOLATION
%
%   Y = MY_INTERPF1M(XX,YY,A)
%
%   XX    = Axis vector for the value matrix YY
%   YY    = Value matrix
%   A     = Set of XX-values to be interpolated
%
%   Y     = Set of interpolated values (same structure as A)
%
%   Assumes equally spaced XX
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström


XX1   = repmat(xx1',length(xx2),1);
XLIMl = repmat(xlim(1,:)',1,length(xx1));
XLIMu = repmat(xlim(2,:)',1,length(xx1));

% find grid point just inside lower boundary
[r c] = find(XX1>XLIMl);
[ru in]=unique(r,'first');
Iinl = c(in);
% Iinl = findl(XX1>XLIMl,'first');

% find grid point just inside upper boundary
[r c] = find(XX1<XLIMu);
[ru in]=unique(r,'last');
Iinu = c(in);
% Iinu = findl(XX1<XLIMu,'last');

% Iinl = find(repmat(xx1',7,1)>repmat(xlim(:,1),1,31),1,'first');
% Iinu = find(xx<xlim(2),1,'last');
xliml = xlim(1,:);
xlimu = xlim(2,:);
yliml = ylim(1,:);
ylimu = ylim(2,:);
xx1l  = xx1(Iinl);
xx1u  = xx1(Iinu);

% find interpolation points between lower boundary and closest grid point
Ibel = find(A1 >=  dpm_interpn(xx2,xliml,A2) & A1 <  dpm_interpn(xx2,xx1l,A2));
% find interpolation points between upper boundary and closest grid point
Ibeu = find(A1 <=  dpm_interpn(xx2,xlimu,A2) & A1 >  dpm_interpn(xx2,xx1u,A2));
% find interpolation points outside boundary
Iout = find(A1 <   dpm_interpn(xx2,xliml,A2) | A1 >  dpm_interpn(xx2,xlimu,A2));
% find interpolation points inside boundary
Iin  = find(A1 >=  dpm_interpn(xx2,xx1l,A2)  & A1 <= dpm_interpn(xx2,xx1u,A2));
% find interpolation points between lower and upper boundary
Iinx = find(A1 >=  dpm_interpn(xx2,xliml,A2) & A1 <= dpm_interpn(xx2,xlimu,A2));

% initialize output
y = nan(size(A1));

% interpolate as usual with interior points
if ~isempty(Iin)
    y(Iin) = dpm_interpn(xx2,xx1,YY,A2(Iin),A1(Iin));
end
% set outside points to inf
if ~isempty(Iout)
    y(Iout)= myInf;
end
% if there are grid points between boundaries
% if ~isempty(find(xx1<xlimu & xx1>xliml,1))
if ~isempty(min(Iinu - Iinl)>0)
    % interpolate points between lower boundary and closest feasible grid point
    if ~isempty(Ibel)
        Xl = dpm_interpn(xx2,xliml,A2(Ibel));
        Xu = dpm_interpn(xx2,xx1l,A2(Ibel));
        Yl = dpm_interpn(xx2,yliml,Xl);
        Yu = dpm_interpn(xx2,YY(dpm_sub2ind(size(YY),Iinl,(1:length(xx2))')),Xu);
        y(Ibel) = Yl + (A1(Ibel) - Xl)./(Xu-Xl).*(Yu-Yl);
        %         y(Ibel) = dpm_interpn([xliml xx1l],[yliml yy(Iinl)],A(Ibel));
    end
    % interpolate points between upper boundary and closest feasible grid point
    if ~isempty(Ibeu)
        Xu = dpm_interpn(xx2,xlimu,A2(Ibeu));
        Xl = dpm_interpn(xx2,xx1u,A2(Ibeu));
        Yu = dpm_interpn(xx2,ylimu,Xu);
        Yl = dpm_interpn(xx2,YY(dpm_sub2ind(size(YY),Iinu,(1:length(xx2))')),Xl);
        y(Ibeu) = Yl + (A1(Ibeu) - Xl)./(Xu-Xl).*(Yu-Yl);
        %         y(Ibeu) = dpm_interpn([xx1u xlimu'],[yy(Iinu) ylim(2)],A2(Ibeu),A1(Ibeu));
    end
else
    % if there are no grid points between boundaries
    y(Iinx) = dpm_interpn(xlim,ylim,A(Iinx));
end
function y         	= dpm_interpf1sbh(xx,yy,a,lim)
%MY_INTERPF1M Computes the 1D interpolation for the given set A
%   using the function YY(xx). USES EXTRAPOLATION
%
%   Y = MY_INTERPF1M(XX,YY,A)
%
%   XX    = Axis vector for the value matrix YY
%   YY    = Value matrix
%   A     = Set of XX-values to be interpolated
%
%   Y     = Set of interpolated values (same structure as A)
%
%   Assumes equally spaced XX
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström
xlu = find(xx >  lim(1),1,'first');
xll = find(xx <= lim(1),1,'last');
xuu = find(xx >= lim(2),1,'first');
xul = find(xx <  lim(2),1,'last');

% if a is between lower limit and regular grid
if a <= lim(1)
    y = yy(xll);
    % if a is outside upper limit
elseif a >= lim(2)
    y = yy(xuu);
    % if a is inside limits and within regular grid
elseif a < xx(xlu) && a > lim(1)
%     % if close to engine off
%     if yy(xlu) == 1 || yy(xll) == 1
%         [tmp ind] = min(abs([xx(xlu) lim(1)]-a));
%         ytmp = [yy(xlu) yy(xll)];
%         y = ytmp(ind);
%     else
        dy  = yy(xlu)-yy(xll);
        dx  = xx(xlu)-lim(1);
        y     = (a-lim(1))*dy/dx + yy(xll);
%     end
% if a is between upper limit and regular grid
elseif  a < lim(2) && a > xx(xul)
%     % if close to engine off
%     if yy(xuu) == 1 || yy(xul) == 1
%         [tmp ind] = min(abs([lim(2) xx(xul)]-a));
%         ytmp = [yy(xuu) yy(xul)];
%         y = ytmp(ind);
%     else
        dy  = yy(xuu)-yy(xul);
        dx  = lim(2)-xx(xul);
        y     = (a-xx(xul))*dy/dx + yy(xul);
%     end
    % if a is outside lower limit
else
%     % if close to engine off
%     il = find(xx<a,1,'last');
%     iu = find(xx>a,1,'first');
%     if yy(il) == 1 || yy(iu)==1
%         ix     = round(dpm_interpn([xx(il) xx(iu)],[1 2],a));
%         yy = [yy(il) yy(iu)];
%         y = yy(ix);
%     else
        y     = dpm_interpn(xx,yy,a);
%     end
end
function y         	= dpm_interpf2sbh(xx1,xx2,YY,a1,a2,lim)
%MY_INTERPF1M Computes the 1D interpolation for the given set A
%   using the function YY(xx). USES EXTRAPOLATION
%
%   Y = MY_INTERPF1M(XX,YY,A)
%
%   XX    = Axis vector for the value matrix YY
%   YY    = Value matrix
%   A     = Set of XX-values to be interpolated
%
%   Y     = Set of interpolated values (same structure as A)
%
%   Assumes equally spaced XX
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström
lim2(1) = dpm_interpn(xx2,lim(1,:),a2);
lim2(2) = dpm_interpn(xx2,lim(2,:),a2);
yy      = dpm_interpn(xx2,xx1,YY,a2.*ones(size(xx1)),xx1);
y       = dpm_interpf1sbh(xx1,yy,a1,lim2);
function y         	= dpm_interpn(varargin)
switch((nargin-1)/2)
    case 1
        xx1 = varargin{1};
        YY  = varargin{2};
        A1  = varargin{3};

        Ars = [reshape(A1, [numel(A1) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);

        h = max([max(diff(xx{1}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);

        yy(:,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1)));
        yy(:,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);

        yi = yy(:,1,1) + da(:,1).*(yy(:,2,1) - yy(:,1,1));

        y = reshape(yi,size(A1));

    case 2
        xx2 = varargin{1};
        xx1 = varargin{2};
        YY  = varargin{3};
        A2  = varargin{4};
        A1  = varargin{5};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);

        yy(:,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1)));
        yy(:,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1)));
        yy(:,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2)));
        yy(:,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);

        v1(:,1) = yy(:,1,1) + da(:,1).*(yy(:,2,1) - yy(:,1,1));
        v1(:,2) = yy(:,1,2) + da(:,1).*(yy(:,2,2) - yy(:,1,2));

        yi = da(:,2).*(v1(:,2) - v1(:,1)) + v1(:,1);

        y = reshape(yi,size(A1));
    case 3
        xx3 = varargin{1};
        xx2 = varargin{2};
        xx1 = varargin{3};
        YY  = varargin{4};
        A3  = varargin{5};
        A2  = varargin{6};
        A1  = varargin{7};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);

        yy(:,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1)));
        yy(:,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1)));
        yy(:,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1)));
        yy(:,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1)));
        yy(:,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2)));
        yy(:,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2)));
        yy(:,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2)));
        yy(:,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);

        
        v2(:,1,1) = yy(:,1,1,1) + da(:,3).*(yy(:,1,1,2) - yy(:,1,1,1));
        v2(:,2,1) = yy(:,2,1,1) + da(:,3).*(yy(:,2,1,2) - yy(:,2,1,1));
        v2(:,1,2) = yy(:,1,2,1) + da(:,3).*(yy(:,1,2,2) - yy(:,1,2,1));
        v2(:,2,2) = yy(:,2,2,1) + da(:,3).*(yy(:,2,2,2) - yy(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);        
        
        
        
        
        
        y = reshape(yi,size(A1));
    case 4
        xx4 = varargin{1};
        xx3 = varargin{2};
        xx2 = varargin{3};
        xx1 = varargin{4};
        YY  = varargin{5};
        A4  = varargin{6};
        A3  = varargin{7};
        A2  = varargin{8};
        A1  = varargin{9};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);

        yy(:,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1)));
        yy(:,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1)));
        yy(:,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1)));
        yy(:,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1)));
        yy(:,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1)));
        yy(:,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1)));
        yy(:,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1)));
        yy(:,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1)));
        yy(:,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2)));
        yy(:,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2)));
        yy(:,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2)));
        yy(:,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2)));
        yy(:,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2)));
        yy(:,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2)));
        yy(:,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2)));
        yy(:,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);

        
        v3(:,1,1,1) = yy(:,1,1,1,1) + da(:,4).*(yy(:,1,1,1,2) - yy(:,1,1,1,1));
        v3(:,2,1,1) = yy(:,2,1,1,1) + da(:,4).*(yy(:,2,1,1,2) - yy(:,2,1,1,1));
        v3(:,1,2,1) = yy(:,1,2,1,1) + da(:,4).*(yy(:,1,2,1,2) - yy(:,1,2,1,1));
        v3(:,2,2,1) = yy(:,2,2,1,1) + da(:,4).*(yy(:,2,2,1,2) - yy(:,2,2,1,1));
        v3(:,1,1,2) = yy(:,1,1,2,1) + da(:,4).*(yy(:,1,1,2,2) - yy(:,1,1,2,1));
        v3(:,2,1,2) = yy(:,2,1,2,1) + da(:,4).*(yy(:,2,1,2,2) - yy(:,2,1,2,1));
        v3(:,1,2,2) = yy(:,1,2,2,1) + da(:,4).*(yy(:,1,2,2,2) - yy(:,1,2,2,1));
        v3(:,2,2,2) = yy(:,2,2,2,1) + da(:,4).*(yy(:,2,2,2,2) - yy(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);          
        
        
        
        y = reshape(yi,size(A1));
    case 5
        xx5 = varargin{1};
        xx4 = varargin{2};
        xx3 = varargin{3};
        xx2 = varargin{4};
        xx1 = varargin{5};
        YY  = varargin{6};
        A5  = varargin{7};
        A4  = varargin{8};
        A3  = varargin{9};
        A2  = varargin{10};
        A1  = varargin{11};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1]) reshape(A5, [numel(A5) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);
        xx{5}  = reshape(xx5,[numel(xx5),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4})) max(diff(xx{5}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});
        Ars(Ars(:,5)<min(xx{5}),5) = min(xx{5});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});
        Ars(Ars(:,5)>max(xx{5}),5) = max(xx{5});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,1)  = 1 + floor(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,2) = 1 +  ceil(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);

        yy(:,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1)));
        yy(:,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1)));
        yy(:,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1)));
        yy(:,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1)));
        yy(:,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2)));
        yy(:,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2)));
        yy(:,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2)));
        yy(:,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2)));
        yy(:,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2)));

        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);
        da(:,5) = (Ars(:,5)-xx{5}(ind(:,5,1)))/h(5);

        v4(:,1,1,1,1) = yy(:,1,1,1,1,1) + da(:,5).*(yy(:,1,1,1,1,2) - yy(:,1,1,1,1,1));
        v4(:,2,1,1,1) = yy(:,2,1,1,1,1) + da(:,5).*(yy(:,2,1,1,1,2) - yy(:,2,1,1,1,1));
        v4(:,1,2,1,1) = yy(:,1,2,1,1,1) + da(:,5).*(yy(:,1,2,1,1,2) - yy(:,1,2,1,1,1));
        v4(:,2,2,1,1) = yy(:,2,2,1,1,1) + da(:,5).*(yy(:,2,2,1,1,2) - yy(:,2,2,1,1,1));
        v4(:,1,1,2,1) = yy(:,1,1,2,1,1) + da(:,5).*(yy(:,1,1,2,1,2) - yy(:,1,1,2,1,1));
        v4(:,2,1,2,1) = yy(:,2,1,2,1,1) + da(:,5).*(yy(:,2,1,2,1,2) - yy(:,2,1,2,1,1));
        v4(:,1,2,2,1) = yy(:,1,2,2,1,1) + da(:,5).*(yy(:,1,2,2,1,2) - yy(:,1,2,2,1,1));
        v4(:,2,2,2,1) = yy(:,2,2,2,1,1) + da(:,5).*(yy(:,2,2,2,1,2) - yy(:,2,2,2,1,1));
        v4(:,1,1,1,2) = yy(:,1,1,1,2,1) + da(:,5).*(yy(:,1,1,1,2,2) - yy(:,1,1,1,2,1));
        v4(:,2,1,1,2) = yy(:,2,1,1,2,1) + da(:,5).*(yy(:,2,1,1,2,2) - yy(:,2,1,1,2,1));
        v4(:,1,2,1,2) = yy(:,1,2,1,2,1) + da(:,5).*(yy(:,1,2,1,2,2) - yy(:,1,2,1,2,1));
        v4(:,2,2,1,2) = yy(:,2,2,1,2,1) + da(:,5).*(yy(:,2,2,1,2,2) - yy(:,2,2,1,2,1));
        v4(:,1,1,2,2) = yy(:,1,1,2,2,1) + da(:,5).*(yy(:,1,1,2,2,2) - yy(:,1,1,2,2,1));
        v4(:,2,1,2,2) = yy(:,2,1,2,2,1) + da(:,5).*(yy(:,2,1,2,2,2) - yy(:,2,1,2,2,1));
        v4(:,1,2,2,2) = yy(:,1,2,2,2,1) + da(:,5).*(yy(:,1,2,2,2,2) - yy(:,1,2,2,2,1));
        v4(:,2,2,2,2) = yy(:,2,2,2,2,1) + da(:,5).*(yy(:,2,2,2,2,2) - yy(:,2,2,2,2,1));


        v3(:,1,1,1) = v4(:,1,1,1,1) + da(:,4).*(v4(:,1,1,1,2) - v4(:,1,1,1,1));
        v3(:,2,1,1) = v4(:,2,1,1,1) + da(:,4).*(v4(:,2,1,1,2) - v4(:,2,1,1,1));
        v3(:,1,2,1) = v4(:,1,2,1,1) + da(:,4).*(v4(:,1,2,1,2) - v4(:,1,2,1,1));
        v3(:,2,2,1) = v4(:,2,2,1,1) + da(:,4).*(v4(:,2,2,1,2) - v4(:,2,2,1,1));
        v3(:,1,1,2) = v4(:,1,1,2,1) + da(:,4).*(v4(:,1,1,2,2) - v4(:,1,1,2,1));
        v3(:,2,1,2) = v4(:,2,1,2,1) + da(:,4).*(v4(:,2,1,2,2) - v4(:,2,1,2,1));
        v3(:,1,2,2) = v4(:,1,2,2,1) + da(:,4).*(v4(:,1,2,2,2) - v4(:,1,2,2,1));
        v3(:,2,2,2) = v4(:,2,2,2,1) + da(:,4).*(v4(:,2,2,2,2) - v4(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);            
        
        y = reshape(yi,size(A1));
    case 6
        xx6 = varargin{1};
        xx5 = varargin{2};
        xx4 = varargin{3};
        xx3 = varargin{4};
        xx2 = varargin{5};
        xx1 = varargin{6};
        YY  = varargin{7};
        A6  = varargin{8};
        A5  = varargin{9};
        A4  = varargin{10};
        A3  = varargin{11};
        A2  = varargin{12};
        A1  = varargin{13};

        Ars = [reshape(A1, [numel(A1) ,1]) reshape(A2, [numel(A2) ,1]) reshape(A3, [numel(A3) ,1]) reshape(A4, [numel(A4) ,1]) reshape(A5, [numel(A5) ,1]) reshape(A6, [numel(A6) ,1])];
        xx{1}  = reshape(xx1,[numel(xx1),1]);
        xx{2}  = reshape(xx2,[numel(xx2),1]);
        xx{3}  = reshape(xx3,[numel(xx3),1]);
        xx{4}  = reshape(xx4,[numel(xx4),1]);
        xx{5}  = reshape(xx5,[numel(xx5),1]);
        xx{6}  = reshape(xx6,[numel(xx6),1]);

        h = max([max(diff(xx{1})) max(diff(xx{2})) max(diff(xx{3})) max(diff(xx{4})) max(diff(xx{5})) max(diff(xx{6}))],eps);

        Ars(Ars(:,1)<min(xx{1}),1) = min(xx{1});
        Ars(Ars(:,2)<min(xx{2}),2) = min(xx{2});
        Ars(Ars(:,3)<min(xx{3}),3) = min(xx{3});
        Ars(Ars(:,4)<min(xx{4}),4) = min(xx{4});
        Ars(Ars(:,5)<min(xx{5}),5) = min(xx{5});
        Ars(Ars(:,6)<min(xx{6}),6) = min(xx{6});

        Ars(Ars(:,1)>max(xx{1}),1) = max(xx{1});
        Ars(Ars(:,2)>max(xx{2}),2) = max(xx{2});
        Ars(Ars(:,3)>max(xx{3}),3) = max(xx{3});
        Ars(Ars(:,4)>max(xx{4}),4) = max(xx{4});
        Ars(Ars(:,5)>max(xx{5}),5) = max(xx{5});
        Ars(Ars(:,6)>max(xx{6}),6) = max(xx{6});

        ind(:,1,1)  = 1 + floor(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,1)  = 1 + floor(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,1)  = 1 + floor(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,1)  = 1 + floor(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,1)  = 1 + floor(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);
        ind(:,6,1)  = 1 + floor(round((Ars(:,6)-xx{6}(1))/h(6)*1e8)*1e-8);

        ind(:,1,2) = 1 +  ceil(round((Ars(:,1)-xx{1}(1))/h(1)*1e8)*1e-8);
        ind(:,2,2) = 1 +  ceil(round((Ars(:,2)-xx{2}(1))/h(2)*1e8)*1e-8);
        ind(:,3,2) = 1 +  ceil(round((Ars(:,3)-xx{3}(1))/h(3)*1e8)*1e-8);
        ind(:,4,2) = 1 +  ceil(round((Ars(:,4)-xx{4}(1))/h(4)*1e8)*1e-8);
        ind(:,5,2) = 1 +  ceil(round((Ars(:,5)-xx{5}(1))/h(5)*1e8)*1e-8);
        ind(:,6,2) = 1 +  ceil(round((Ars(:,6)-xx{6}(1))/h(6)*1e8)*1e-8);

        yy(:,1,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,1,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,2,1,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,1,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,1,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,2,2,2,2,1,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,1)));
        yy(:,1,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,1,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,2,1,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,1,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,1,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,2,2,2,2,2,1)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,1)));
        yy(:,1,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,1,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,2,1,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,1,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,1,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,2,2,2,2,1,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,1),ind(:,6,2)));
        yy(:,1,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,1,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,2,1,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,1),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,1,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,1),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,1,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,1),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,1,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,1),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        yy(:,2,2,2,2,2,2)    = YY(dpm_sub2ind(size(YY),ind(:,1,2),ind(:,2,2),ind(:,3,2),ind(:,4,2),ind(:,5,2),ind(:,6,2)));
        
        da(:,1) = (Ars(:,1)-xx{1}(ind(:,1,1)))/h(1);
        da(:,2) = (Ars(:,2)-xx{2}(ind(:,2,1)))/h(2);
        da(:,3) = (Ars(:,3)-xx{3}(ind(:,3,1)))/h(3);
        da(:,4) = (Ars(:,4)-xx{4}(ind(:,4,1)))/h(4);
        da(:,5) = (Ars(:,5)-xx{5}(ind(:,5,1)))/h(5);
        da(:,6) = (Ars(:,6)-xx{6}(ind(:,6,1)))/h(6);

        v5(:,1,1,1,1,1) = yy(:,1,1,1,1,1,1) + da(:,6).*(yy(:,1,1,1,1,1,2) - yy(:,1,1,1,1,1,1));
        v5(:,2,1,1,1,1) = yy(:,2,1,1,1,1,1) + da(:,6).*(yy(:,2,1,1,1,1,2) - yy(:,2,1,1,1,1,1));
        v5(:,1,2,1,1,1) = yy(:,1,2,1,1,1,1) + da(:,6).*(yy(:,1,2,1,1,1,2) - yy(:,1,2,1,1,1,1));
        v5(:,2,2,1,1,1) = yy(:,2,2,1,1,1,1) + da(:,6).*(yy(:,2,2,1,1,1,2) - yy(:,2,2,1,1,1,1));
        v5(:,1,1,2,1,1) = yy(:,1,1,2,1,1,1) + da(:,6).*(yy(:,1,1,2,1,1,2) - yy(:,1,1,2,1,1,1));
        v5(:,2,1,2,1,1) = yy(:,2,1,2,1,1,1) + da(:,6).*(yy(:,2,1,2,1,1,2) - yy(:,2,1,2,1,1,1));
        v5(:,1,2,2,1,1) = yy(:,1,2,2,1,1,1) + da(:,6).*(yy(:,1,2,2,1,1,2) - yy(:,1,2,2,1,1,1));
        v5(:,2,2,2,1,1) = yy(:,2,2,2,1,1,1) + da(:,6).*(yy(:,2,2,2,1,1,2) - yy(:,2,2,2,1,1,1));
        v5(:,1,1,1,2,1) = yy(:,1,1,1,2,1,1) + da(:,6).*(yy(:,1,1,1,2,1,2) - yy(:,1,1,1,2,1,1));
        v5(:,2,1,1,2,1) = yy(:,2,1,1,2,1,1) + da(:,6).*(yy(:,2,1,1,2,1,2) - yy(:,2,1,1,2,1,1));
        v5(:,1,2,1,2,1) = yy(:,1,2,1,2,1,1) + da(:,6).*(yy(:,1,2,1,2,1,2) - yy(:,1,2,1,2,1,1));
        v5(:,2,2,1,2,1) = yy(:,2,2,1,2,1,1) + da(:,6).*(yy(:,2,2,1,2,1,2) - yy(:,2,2,1,2,1,1));
        v5(:,1,1,2,2,1) = yy(:,1,1,2,2,1,1) + da(:,6).*(yy(:,1,1,2,2,1,2) - yy(:,1,1,2,2,1,1));
        v5(:,2,1,2,2,1) = yy(:,2,1,2,2,1,1) + da(:,6).*(yy(:,2,1,2,2,1,2) - yy(:,2,1,2,2,1,1));
        v5(:,1,2,2,2,1) = yy(:,1,2,2,2,1,1) + da(:,6).*(yy(:,1,2,2,2,1,2) - yy(:,1,2,2,2,1,1));
        v5(:,2,2,2,2,1) = yy(:,2,2,2,2,1,1) + da(:,6).*(yy(:,2,2,2,2,1,2) - yy(:,2,2,2,2,1,1));
        v5(:,1,1,1,1,2) = yy(:,1,1,1,1,2,1) + da(:,6).*(yy(:,1,1,1,1,2,2) - yy(:,1,1,1,1,2,1));
        v5(:,2,1,1,1,2) = yy(:,2,1,1,1,2,1) + da(:,6).*(yy(:,2,1,1,1,2,2) - yy(:,2,1,1,1,2,1));
        v5(:,1,2,1,1,2) = yy(:,1,2,1,1,2,1) + da(:,6).*(yy(:,1,2,1,1,2,2) - yy(:,1,2,1,1,2,1));
        v5(:,2,2,1,1,2) = yy(:,2,2,1,1,2,1) + da(:,6).*(yy(:,2,2,1,1,2,2) - yy(:,2,2,1,1,2,1));
        v5(:,1,1,2,1,2) = yy(:,1,1,2,1,2,1) + da(:,6).*(yy(:,1,1,2,1,2,2) - yy(:,1,1,2,1,2,1));
        v5(:,2,1,2,1,2) = yy(:,2,1,2,1,2,1) + da(:,6).*(yy(:,2,1,2,1,2,2) - yy(:,2,1,2,1,2,1));
        v5(:,1,2,2,1,2) = yy(:,1,2,2,1,2,1) + da(:,6).*(yy(:,1,2,2,1,2,2) - yy(:,1,2,2,1,2,1));
        v5(:,2,2,2,1,2) = yy(:,2,2,2,1,2,1) + da(:,6).*(yy(:,2,2,2,1,2,2) - yy(:,2,2,2,1,2,1));
        v5(:,1,1,1,2,2) = yy(:,1,1,1,2,2,1) + da(:,6).*(yy(:,1,1,1,2,2,2) - yy(:,1,1,1,2,2,1));
        v5(:,2,1,1,2,2) = yy(:,2,1,1,2,2,1) + da(:,6).*(yy(:,2,1,1,2,2,2) - yy(:,2,1,1,2,2,1));
        v5(:,1,2,1,2,2) = yy(:,1,2,1,2,2,1) + da(:,6).*(yy(:,1,2,1,2,2,2) - yy(:,1,2,1,2,2,1));
        v5(:,2,2,1,2,2) = yy(:,2,2,1,2,2,1) + da(:,6).*(yy(:,2,2,1,2,2,2) - yy(:,2,2,1,2,2,1));
        v5(:,1,1,2,2,2) = yy(:,1,1,2,2,2,1) + da(:,6).*(yy(:,1,1,2,2,2,2) - yy(:,1,1,2,2,2,1));
        v5(:,2,1,2,2,2) = yy(:,2,1,2,2,2,1) + da(:,6).*(yy(:,2,1,2,2,2,2) - yy(:,2,1,2,2,2,1));
        v5(:,1,2,2,2,2) = yy(:,1,2,2,2,2,1) + da(:,6).*(yy(:,1,2,2,2,2,2) - yy(:,1,2,2,2,2,1));
        v5(:,2,2,2,2,2) = yy(:,2,2,2,2,2,1) + da(:,6).*(yy(:,2,2,2,2,2,2) - yy(:,2,2,2,2,2,1));

        v4(:,1,1,1,1) = v5(:,1,1,1,1,1) + da(:,5).*(v5(:,1,1,1,1,2) - v5(:,1,1,1,1,1));
        v4(:,2,1,1,1) = v5(:,2,1,1,1,1) + da(:,5).*(v5(:,2,1,1,1,2) - v5(:,2,1,1,1,1));
        v4(:,1,2,1,1) = v5(:,1,2,1,1,1) + da(:,5).*(v5(:,1,2,1,1,2) - v5(:,1,2,1,1,1));
        v4(:,2,2,1,1) = v5(:,2,2,1,1,1) + da(:,5).*(v5(:,2,2,1,1,2) - v5(:,2,2,1,1,1));
        v4(:,1,1,2,1) = v5(:,1,1,2,1,1) + da(:,5).*(v5(:,1,1,2,1,2) - v5(:,1,1,2,1,1));
        v4(:,2,1,2,1) = v5(:,2,1,2,1,1) + da(:,5).*(v5(:,2,1,2,1,2) - v5(:,2,1,2,1,1));
        v4(:,1,2,2,1) = v5(:,1,2,2,1,1) + da(:,5).*(v5(:,1,2,2,1,2) - v5(:,1,2,2,1,1));
        v4(:,2,2,2,1) = v5(:,2,2,2,1,1) + da(:,5).*(v5(:,2,2,2,1,2) - v5(:,2,2,2,1,1));
        v4(:,1,1,1,2) = v5(:,1,1,1,2,1) + da(:,5).*(v5(:,1,1,1,2,2) - v5(:,1,1,1,2,1));
        v4(:,2,1,1,2) = v5(:,2,1,1,2,1) + da(:,5).*(v5(:,2,1,1,2,2) - v5(:,2,1,1,2,1));
        v4(:,1,2,1,2) = v5(:,1,2,1,2,1) + da(:,5).*(v5(:,1,2,1,2,2) - v5(:,1,2,1,2,1));
        v4(:,2,2,1,2) = v5(:,2,2,1,2,1) + da(:,5).*(v5(:,2,2,1,2,2) - v5(:,2,2,1,2,1));
        v4(:,1,1,2,2) = v5(:,1,1,2,2,1) + da(:,5).*(v5(:,1,1,2,2,2) - v5(:,1,1,2,2,1));
        v4(:,2,1,2,2) = v5(:,2,1,2,2,1) + da(:,5).*(v5(:,2,1,2,2,2) - v5(:,2,1,2,2,1));
        v4(:,1,2,2,2) = v5(:,1,2,2,2,1) + da(:,5).*(v5(:,1,2,2,2,2) - v5(:,1,2,2,2,1));
        v4(:,2,2,2,2) = v5(:,2,2,2,2,1) + da(:,5).*(v5(:,2,2,2,2,2) - v5(:,2,2,2,2,1));

        
        v3(:,1,1,1) = v4(:,1,1,1,1) + da(:,4).*(v4(:,1,1,1,2) - v4(:,1,1,1,1));
        v3(:,2,1,1) = v4(:,2,1,1,1) + da(:,4).*(v4(:,2,1,1,2) - v4(:,2,1,1,1));
        v3(:,1,2,1) = v4(:,1,2,1,1) + da(:,4).*(v4(:,1,2,1,2) - v4(:,1,2,1,1));
        v3(:,2,2,1) = v4(:,2,2,1,1) + da(:,4).*(v4(:,2,2,1,2) - v4(:,2,2,1,1));
        v3(:,1,1,2) = v4(:,1,1,2,1) + da(:,4).*(v4(:,1,1,2,2) - v4(:,1,1,2,1));
        v3(:,2,1,2) = v4(:,2,1,2,1) + da(:,4).*(v4(:,2,1,2,2) - v4(:,2,1,2,1));
        v3(:,1,2,2) = v4(:,1,2,2,1) + da(:,4).*(v4(:,1,2,2,2) - v4(:,1,2,2,1));
        v3(:,2,2,2) = v4(:,2,2,2,1) + da(:,4).*(v4(:,2,2,2,2) - v4(:,2,2,2,1));
        
        v2(:,1,1) = v3(:,1,1,1) + da(:,3).*(v3(:,1,1,2) - v3(:,1,1,1));
        v2(:,2,1) = v3(:,2,1,1) + da(:,3).*(v3(:,2,1,2) - v3(:,2,1,1));
        v2(:,1,2) = v3(:,1,2,1) + da(:,3).*(v3(:,1,2,2) - v3(:,1,2,1));
        v2(:,2,2) = v3(:,2,2,1) + da(:,3).*(v3(:,2,2,2) - v3(:,2,2,1)); 
        
        v1(:,1) = v2(:,1,1) + da(:,2).*(v2(:,1,2) - v2(:,1,1));
        v1(:,2) = v2(:,2,1) + da(:,2).*(v2(:,2,2) - v2(:,2,1));

        yi = da(:,1).*(v1(:,2) - v1(:,1)) + v1(:,1);            
        
        y = reshape(yi,size(A1));
    otherwise
        error('DPM:Internal','Too many states or inputs: contact the author of DPM')
end
function S        	= dpm_mergestruct(S1,S2)
% WARNING: might cause memory problems!
try
    S = S1;
    names = fieldnames(S1);
    for i=1:length(names)
        if isstruct(S1.(names{i}))
            S.(names{i}) = dpm_mergestruct(S1.(names{i}), S2.(names{i}));
        elseif iscell(S1.(names{i}))
            for j=1:numel(S1.(names{i}))
                S.(names{i}){j} = [S1.(names{i}){j} S2.(names{i}){j}];
            end
        elseif isnumeric(S1.(names{i})) || islogical(S1.(names{i}))
            S.(names{i}) = [S1.(names{i}) S2.(names{i})];
        else
            S.(names{i}) = S2.(names{i});
        end
    end
catch
    error('mergestruct: S1 and S2 have different structures.')
end
function inp       	= dpm_get_empty_inp(grd,dis,options)
if strcmp(options,'nan')
    value = nan;
elseif strcmp(options,'zero')
    value = 0;
elseif strcmp(options,'inf')
    value = inf;
end

for i=1:length(grd.Nx)
    inp.X{i} = value;
end
for i=1:length(grd.Nu)
    inp.U{i} = value;
end
for i=1:length(dis.W)
    inp.W{i} = value;
end
inp.Ts = dis.Ts;
function out       	= dpm_get_empty_out(model,inp,par,grd,options)
%GET_EMPTY_OUT   Gets an empty result struct
%   OUT   = GET_EMPTY_OUT(MODEL,OPTIONS)
%   Gets an empty output struct from model.
%
%   See also dpm_forward_sim
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström
if ~exist('options')
    options = 'nan';
end

[X C I out] = feval(model,inp,par);
out.X = X;
out.C = C;
out.I = I;
out = dpm_setallfield(out,options);

function S         	= dpm_setallfield(S1,options)

try
    if strcmp(options,'nan')
        value = nan;
    elseif strcmp(options,'zero')
        value = 0;
    elseif strcmp(options,'inf')
        value = inf;
    end

    S = S1;
    names = fieldnames(S1);
    for i=1:length(names)
        if isstruct(S1.(names{i}))
            S.(names{i}) = dpm_setallfield(S1.(names{i}), options);
        elseif iscell(S1.(names{i}))
            for j=1:numel(S1.(names{i}))
                S.(names{i}){j} = value;
            end
        elseif isnumeric(S1.(names{i})) || islogical(S1.(names{i}))
            S.(names{i}) = value;
        else
            S.(names{i}) = S1.(names{i});
        end
    end
catch
    error('mergestruct: S1 and S2 have different structures.')
end

function [X C I out]= dpm_model_inv(inp,par)
inpo = inp;
iterations = 0;
dSOC = inf;

while max(abs(reshape(dSOC,1,numel(dSOC)))) > par.options.Tol && iterations < par.options.Iter
%     [out.X out.C out.I] = eval([par.model '(inp,par);']);
    [X C I] = feval(par.model,inp,par);
    dSOC   = X{1} - inpo.X{1};
    inp.X{1} = inp.X{1} - dSOC;
    iterations = iterations+1;
end
% out.I = bitor(out.I,bitor(out.X{1}<lim.Xn{1}.lo, out.X{1}>lim.Xn{1}.hi));
X = inp.X;
C{2} = C{1};
C{1} = (X{1}-inpo.X{1});
out = [];

function str       	= dpm_code(s,v,m)
% v = [1 2 3];
% s = 'inp.U{*}';
% str = '';
if ~exist('m','var')
    m='';
end
v    = reshape(v,1,numel(v));
str  = repmat([s m],1,length(v));
str  = str(1:end-length(m));
star = strfind(str, '#');
str(star) = strrep(int2str(v), ' ', '');

function in         = dpm_findu(A,vec)
da = A(2)-A(1);
in = 1+ceil((vec-A(1))./da);
in = max(in,1);
in = min(in,length(A));
function in         = dpm_findl(A,vec)
da = A(2)-A(1);
in = 1+floor((vec-A(1))./da);
in = max(in,1);
in = min(in,length(A));
function [ind col]  = dpm_sub2indr(sze,vl,vu,dim)

ind0 = [0:sze(dim)-1].*sze(1);
ind = reshape([(vl+ind0); (vu+ind0)],1,numel(ind0)*2);
indstr = num2str(ind');
ind = eval(['[' reshape([indstr repmat(': ',1,length(vl))']',1,numel([indstr repmat(': ',1,length(vl))']')) ']']);
col = ceil((ind)/sze(1));
function ndx        = dpm_sub2ind(siz,varargin)

siz = double(siz);

if length(siz) ~= nargin-1
    %Adjust input
    if length(siz)<nargin-1
        %Adjust for trailing singleton dimensions
        siz = [siz ones(1,nargin-length(siz)-1)];
    else
        %Adjust for linear indexing on last element
        siz = [siz(1:nargin-2) prod(siz(nargin-1:end))];
    end
end

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:length(siz),
    v = varargin{i};
    ndx = ndx + (v-1)*k(i);
end
function c          = dpm_sizecmp(a,b)
sa = size(a);
sb = size(b);
c = numel(sa)==numel(sb) & numel(a) == numel(b) & sum(sa==sb);

function sze        = get_size(current_grd)
sze = [];
for i=1:length(current_grd.X)
    sze = [sze length(current_grd.X{i})];
end
if length(sze)==1
    sze = [sze 1];
end

function grd        = input_check_grd(grd,T)

% grd = 
%     Xn: {[1x1 struct]}
%     XN: {[1x1 struct]}
%     X0: {[250]}
%     Nx: {[1x1001 double]}
%     Nu: {[1x1001 double]}
%     Un: {[1x1 struct]}
if T < 1
        error('DPM:Internal','prb.N must be greater than 0')
end
for i=1:length(grd.Nx)
    if grd.Nx{i} < 1
        error('DPM:Internal','grd.Nx{.} must be equal or greater than 1')
    end
    if length(grd.Xn{i}.lo)==1
        grd.Xn{i}.lo = repmat(grd.Xn{i}.lo,1,T+1);
    elseif length(grd.Xn{i}.lo)~=T+1
        error('DPM:Internal','grd.Xn{.}.lo must be a scalar OR have the same length as the problem')
    end
    if length(grd.Xn{i}.hi)==1
        grd.Xn{i}.hi = repmat(grd.Xn{i}.hi,1,T+1);
    elseif length(grd.Xn{i}.hi)~=T+1
        error('DPM:Internal','grd.Xn{.}.hi must be a scalar OR have the same length as the problem')
    end
    if length(grd.Nx{i})==1
        grd.Nx{i} = repmat(grd.Nx{i},1,T+1);
    elseif length(grd.Nx{i})~=T+1
        error('DPM:Internal','grd.Nx{.} must be a scalar OR have the same length as the problem')
    end
end
for i=1:length(grd.Nu)
    if grd.Nu{i} < 1
        error('DPM:Internal','grd.Nu{.} must be equal or greater than 1')
    end
    if length(grd.Un{i}.lo)==1
        grd.Un{i}.lo = repmat(grd.Un{i}.lo,1,T);
    elseif length(grd.Un{i}.lo)~=T
        error('DPM:Internal','grd.Un{.}.lo must be a scalar OR have the same length as the problem')
    end
    if length(grd.Un{i}.hi)==1
        grd.Un{i}.hi = repmat(grd.Un{i}.hi,1,T);
    elseif length(grd.Un{i}.hi)~=T
        error('DPM:Internal','grd.Un{.}.hi must be a scalar OR have the same length as the problem')
    end
    if length(grd.Nu{i})==1
        grd.Nu{i} = repmat(grd.Nu{i},1,T);
    elseif length(grd.Nu{i})~=T
        error('DPM:Internal','grd.Nu{.} must be a scalar OR have the same length as the problem')
    end
end

function              notify_user_of_error(err)

    clear_waitbars();

    switch(err.identifier)
        case 'MATLAB:unassignedOutputs'
            fprintf('DPM:Model function error \n \t Make sure all the output arguments are set in the \n\t model function (X, C, I, out).\n')
        case 'MATLAB:class:SetProhibited'
            fprintf('DPM:Model function error \n \t Make sure all the output arguments are set in the \n\t model function (X, C, I, out).\n')
            fprintf('  or \n \t Make sure the each of the output states X{1}, X{2},..., X{Nx} have \n\t the same dimensions as input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
            fprintf('  or \n \t Make sure the each of the outputs C{1} and I have \n\t the same dimensions as the input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
        case 'MATLAB:cellRefFromNonCell'
            fprintf('DPM:Model function error \n \t Make sure all the output C is a 1x1 cell array and \n\t that the output X is a cell array with as many elements \n\t as state variables.\n')
        case 'MATLAB:dimagree'
            fprintf('DPM:Model function error \n \t Make sure the each of the output states X{1}, X{2},..., X{Nx} have \n\t the same dimensions as input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
        case 'DPM:Internal'
            %fprintf(['DPM:Error \n \t' err.message '\n'])
            fprintf('DPM:Error \n \t Check the model function.\n')
            rethrow(err)
        otherwise
%             err.identifier
            fprintf('DPM:Error \n \t Check the model function.\n')
            rethrow(err)
    end

function clear_waitbars()
        
    % clear all waitbars
    set(0,'ShowHiddenHandles','on')
    handles = get(0,'Children');
    for i=1:length(handles)
        if strcmp(get(handles(i),'name'),'DPM:Waitbar')
            delete(handles(i))
        end
    end
    set(0,'ShowHiddenHandles','off')        