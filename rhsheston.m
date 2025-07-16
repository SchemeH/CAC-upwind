function f = rhsheston(t, y, A, b);
% RHSHESTON - Computes the full rhs of the discretised Heston PDE
%������ɢ��Heston PDE������rhs 
%             y_t = A y + b 
%
% DESCRIPTION:
%   This is used for the matlab function ode15s.
%
% PARAMETERS: 
%   t    - is a scalar of current time.�ǵ�ǰʱ��ı�����
%   y    - is a vector of solution values.�ǽ�ֵ������
%   A    - constant matrix.
%   b    - vector of boundary conditions.
%
% RETURNS:
%   f    - the rhs.

f = A*y+b;
