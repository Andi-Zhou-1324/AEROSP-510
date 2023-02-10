%%%%%
% n = number o f g ri d p oi n t s
% nk = n ( number o f s p e c t r a l p oi n t s ; Note t h a t onl y nk /2 e f f e c t i v e
% p oi n t s due t o p e r i o d i c boundary c o n di ti o n )
% dk = s p a ci n g between p oi n t s i n s p e c t r a l space , taken t o be 1 h e r e
%%% Get the mesh i n s p e c t r a l s p ace

[U,V,W] = read128();

n = size(U,1);

L = 2* pi;
dk = 2* pi /L;
nk = n;
for i =1:nk/2+1
    kx ( i ) = (i - 1)*dk;
end
for i = nk /2+2: nk
    kx ( i ) = -(nk+1-i )* dk ;
end
%%%% To compute F o u ri e r t r an s f o rm , u se f f t n
%% Spectrum i s a s p h e r i c a l a v e r a g e
for k = 1 : nk
    for j = 1 : nk
        for i = 1 : nk
            kk = sqrt ( kx ( i )^2 + kx ( j )^2 + kx ( k )^2 ) ;
            ik = 1 + floor(kk/dk + 0.5);
            if (kk >1e-06)
                spect = 0.5*(abs(Uk(i ,j ,k))^2+ abs(Vk(i ,j ,k ))^2+ abs(Wk(i ,j ,k))^2);
                spectrum(ik) = spectrum (ik) + spect;
            end
        end
    end
end