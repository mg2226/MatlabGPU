function noise_var=set_Q_noise_solve(rule, vk, y,virusobjs, p_theta_eta,Na)

WEIGHT_INDEX = 6; % the index of the abscissa weights in the rule data structure
vo = virusobjs; % use a shorter name inside the function
Nzi = size(rule, 1);
A=0;
B=0;
eta=1;
nu=vo{eta}.nu;

% read tilde b
b_fn = 'callti.out.read_by_C'; % file contains \tilde{b} values computed from mathematica
lmax = max(vo{eta}.clnp.il);
tilde_b = rd_b(b_fn, lmax);

parfor (n=1:Nzi)
    Rabc = euler2R(rule(n, 1:3));
    L=setL_nord(Rabc, vo{eta}.clnp.il, vo{eta}.clnp.in, vk, vo{eta}.Htable, vo{eta}.map_unique2lp,tilde_b);
    Lc = L*vo{eta}.cbar; % \mu = L*c;
    y_Lc = bsxfun(@minus, y, Lc);
    wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(:,n);
    A=A+sum(wri);
    B=B+ sum(y_Lc.^2)*wri;
end
A=A*Na*Na;
noise_var=B/A;
