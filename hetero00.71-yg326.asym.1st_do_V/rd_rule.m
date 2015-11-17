function rule=rd_rule(rule_fn)
%function rule=rd_rule(rule_fn)

if (~exist(rule_fn,'file'))
  error(['rd_rule: fopen: ' rule_fn]);
end

rule = load(rule_fn);