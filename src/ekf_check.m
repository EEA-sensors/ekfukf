%EKF_CHECK   Check derivatives for the EKF's dynamic and measurement
%            models using finite differences.
%
% Syntax:
%  [DF0,DF1,DH0,DH1] = EKF_CHECK(F,DF,)
%

function ekf_check(f,df,index_f,h,hf,index_h,varargin)
  % Calculate function value and derivative
  if isstr(f) | strcmp(class(f),'function_handle')
    yf0 = feval(f,varargin{:});
  else
    yf0 = f(varargin{:});
  end
  if isnumeric(df)
    DF0 = df;
  elseif isstr(df) | strcmp(class(df),'function_handle')
    DF0 = feval(df,varargin{:});
  else
    DF0 = df(varargin{:});
  end  
    
  % Calculate function value and derivative
  if isstr(h) | strcmp(class(h),'function_handle')
    yh0 = feval(h,varargin{:});
  else
    yh0 = h(varargin{:});
  end
  if isnumeric(dh)
    DH0 = dh;
  elseif isstr(dh) | strcmp(class(dh),'function_handle')
    DH0 = feval(dh,varargin{:});
  else
    DH0 = dh(varargin{:});
  end  






