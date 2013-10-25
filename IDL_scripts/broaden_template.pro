function broaden_template, template, velScale, pars



vel = pars[0]/velScale     ; in pixels

sigma = pars[1]/velScale   ; in pixels

dx = ceil(abs(vel)+4*sigma) ; Sample the Gaussian and GH at least to vel+4*sigma

x = range(dx,-dx)           ; Evaluate the Gaussian using steps of 1 pixel.

w = (x - vel)/sigma

w2 = w^2

losvd = exp(-0.5d*w2)/(sqrt(2d*!dpi)*sigma) ; Normalized total(Gaussian)=1

;

; (h5,h6) polynomials from Appendix C of Cappellari et al. (2002)

;

npars = n_elements(pars)

if npars gt 2 then begin

    poly = 1d + pars[2]/Sqrt(3d)*(w*(2d*w2-3d)) $       ; H3

              + pars[3]/Sqrt(24d)*(w2*(4d*w2-12d)+3d)   ; H4

    if npars eq 6 then $

        poly = poly + pars[4]/Sqrt(60d)*(w*(w2*(4d*w2-20d)+15d)) $      ; H5

                    + pars[5]/Sqrt(720d)*(w2*(w2*(8d*w2-60d)+90d)-15d)  ; H6

    losvd = losvd*poly

endif

;

s = size(template)

npix = s[1]

if s[0] eq 2 then ntemp = s[2] else ntemp = 1 ; Number of template spectra

c = dblarr(npix,ntemp,/NOZERO) ; This array is used for estimating predictions

for j=0,ntemp-1 do c[*,j] = convol(template[*,j],losvd,/EDGE_TRUNCATE)

return,c
end
