function NPP=vgpm(chl, par, sst, Dirr)

%  Coded by Stephanie Henson
%  output in mgC/m2/day
% Description:      computes daily primary productivity using
%                      the Behrenfeld-Falkowski (BeFa) algorithm.  The BeFa
%                      algorithm estimates productivity using surface chl
%                      (mg m-3), surface irradiance (Einsteins m-2 d-1),
%                      sea surface temperature (C), and day length (hours).
% 		     Pb_opt is modelled as a polynomial function of SST.
%
%    Input Parameters:
%       chl            Chlorophyll_a surface concentration in milligrams
%                      chlorophyl per cubic meter
%       par            Photosynthetically available radiation in Einsteins per
%                      day per square meter
%       sst            Sea surface temperature in degrees Centigrade
%       Dirr           Length day in decimal hours.
%
%  Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
%       from Satellite-Based Chlorophyll Concentration.  Limnology and
%       Oceanography, Volume 42, Number 1
%
%   NOTES: CHL, PAR and SST should be monthly means, structured as lon X lat
%   X time.  Make sure the CHL, PAR and SST datasets all cover the same time
%   period


for n=1:size(chl,3)

  CHL=chl(:,:,n);
  SST=sst(:,:,n);
  PAR=par(:,:,n);

  if CHL<1
    Ctot=38*(CHL.^0.425); % integrated water column chl
  else
    Ctot=40.2*(CHL.^0.507);
  end

  Zeu=568.2*(Ctot.^-0.746); % euphotic zone depth

  if Zeu>102
    Zeu=200*(Ctot.^-0.293);
  end

  if SST<-1   %% the simple (!) estimate of maximum daily PP
    PBopt=1.13;
  elseif SST>28.5
    PBopt=4;
  else
    PBopt=1.2956 + SST.*2.749e-1 + (SST.^2)*6.17e-2 - (SST.^3)*2.05e-2...
      + (SST.^4)*2.462e-3 - (SST.^5)*1.348e-4 + (SST.^6)*3.4132e-6...
      - (SST.^7)*3.27e-8;
  end

  x=sym(n/12);
  xx=frac(x);
  if xx==0
    N=12;
  else
    N=double(xx)*12;
  end

  DLen=Dirr(:,:,N);  % decimal day length

  NPP(:,:,n)=0.66125*(PBopt.*(PAR./(PAR+4.1)).*CHL.*Zeu.*DLen);
end

NPP(find(NPP==0))=NaN;
