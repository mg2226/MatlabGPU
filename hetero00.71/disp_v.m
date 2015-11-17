function disp_v(verbosity,vlevel,x)
%function disp_v(verbosity,vlevel,x)

if verbosity>=vlevel
  disp(x);
end
