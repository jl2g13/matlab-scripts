% Plotting J - light limitation and temperature dependence term
% J is realised growth rate under given light,T conditions (i.e. assuming
% nutrient replete conditions; units mgC/d?)

% J(I,ChlC,T)

% For observational E(saturation) values see Kirk Table 10.1 p349: E-saturation model values (Ek,Esat) look good


Temp = [-2 0 5 10 20]; %double(-2); % water temperature - For testing J sensitivity to T

for i = 1:length(Temp)    
    
% Temperature depedence

% loc_T = double(0); %double(-2); % water temperature
loc_T = Temp(i); % For testing J sensitivity to T
xvpn = double(0.53); % xvpn = Vpn in paper - max nd phytoplankton growth rate at 0°C. = 0.53 (Yool2010 table 1)
xvpd = double(0.50); % xvpd = Vpd in paper - max d phytoplankton growth rate at 0°C. = 0.50 (Yool2010 table 1)

fun_T = 1.066 ^ (loc_T); % Yool2010 eq14
xvpnT = xvpn * fun_T;
xvpdT = xvpd * fun_T;


% Light field and photodependence (Chl/biomass)

alphan = double(15.0); alphad = double(11.25); % Non-diatom and diatom chl-specific initial slope of P-I curve (Yool2010 Table1)
% eps    = double(0.01257); % Chl:N conversion factor (assume constant C/N = 6.625) (Yool2010 Table 1)
ChlC   = 0:(0.05/200):0.05; % = (chl/ Pn) * eps (Yool2010 eq12) calculated in model from Chl and Pn fields by Yool2010 eq12: in situ published range = 24-1250 (OASIS_MATRAI_2013 slides), model range = 50-200 (model enforced min of 20, Yool 2010 Table 1), satellite = 90.

xpar   = 0:1:200; % Matsuoka2009 fig2a = 0-359 W/m2 (sea surface).

faln = alphan * ChlC; % Yool2010 eq13
fald = alphad * ChlC; % Yool2010 eq20




pos = 1;
    for j = 1:length(faln)
    
% verbose (original fortran translated) and brief directly equivalent for all applications so far (11/4/14)         
        
%         % (1) non-diatoms

             % (a) verbose form

% %         fchn1 = (xvpnT.^2) + ((faln(j).^2) .* (xpar.^2)); % Yool2010 eq15 (denominator)
%           fchn1 = (xvpnT * xvpnT) + (faln(j) .* faln(j) .* xpar .* xpar);
% %         fchn  = xvpnT ./ sqrt(fchn1);
%           if fchn1 > 0.001
%             fchn = xvpnT ./ sqrt(fchn1);
%           else
%             fchn = 0;  
%           end
%           
%         % fjln(:,pos,i)  = fchn .* faln(j) .* xpar; % Yool2010 eq 15 (full)
%         temp = fchn .* faln .* xpar; % Yool2010 eq 15 (full)


        
            % (b) brief form
          
fjln(:,pos,i) = (xvpnT ./ sqrt((xvpnT * xvpnT) + (faln(j) .* faln(j) .* xpar .* xpar))) .* faln(j) .* xpar; % Yool2010 eq 15
          

        % diatoms
        
             % (a) verbose form 
        
% %         fchd1 = (xvpdT.^2) + ((fald(j).^2) .* (xpar.^2)); % Yool2010 eq22 (denominator)
%           fchd1 = (xvpdT * xvpdT) + (fald(j) .* fald(j) .* xpar .* xpar); 
% %         fchd  = xvpdT ./ sqrt(fchd1);
%           if fchd1 > 0.001
%             fchd  = xvpdT ./ sqrt(fchd1);
%           else
%             fchd = 0;
%           end
%                     
%          fjld(:,pos,i)  = fchd .* fald(j) .* xpar; % Yool2010 eq22 (full)

             % (b) brief form
        
          fjld(:,pos,i) = ((xvpdT .* fald(j) .* xpar) ./ sqrt((xvpdT * xvpdT) + (fald(j) .* fald(j) .* xpar .* xpar))); % Yool2010 eq22

        pos = pos + 1;
    end

end

% Plotting

% % Subplot grid of temperature dependence on J term 3D plot
% figure(1);
% subplot(2,4,1); pcolor(xpar,ChlC,fjln(:,:,1)); shading flat; title('Non-diatom J term at T = -2°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,2); pcolor(xpar,ChlC,fjln(:,:,2)); shading flat; title('Non-diatom J term at T = 0°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,3); pcolor(xpar,ChlC,fjln(:,:,3)); shading flat; title('Non-diatom J term at T = 5°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,4); pcolor(xpar,ChlC,fjln(:,:,5)); shading flat; title('Non-diatom J term at T = 20°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,5); pcolor(xpar,ChlC,fjld(:,:,1)); shading flat; title('Diatom J term at T = -2°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,6); pcolor(xpar,ChlC,fjld(:,:,2)); shading flat; title('Diatom J term at T = 0°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,7); pcolor(xpar,ChlC,fjld(:,:,3)); shading flat; title('Diatom J term at T = 5°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% subplot(2,4,8); pcolor(xpar,ChlC,fjld(:,:,5)); shading flat; title('Diatom J term at T = 20°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J');
% 
% % %Figure 1b - caxis=0-1.9 scale made the same for all subplots
% % figure(1);
% % subplot(2,4,1); pcolor(xpar,ChlC,fjln(:,:,1)); shading flat; title('Non-diatom J term at T = -2°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,2); pcolor(xpar,ChlC,fjln(:,:,2)); shading flat; title('Non-diatom J term at T = 0°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,3); pcolor(xpar,ChlC,fjln(:,:,3)); shading flat; title('Non-diatom J term at T = 5°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,4); pcolor(xpar,ChlC,fjln(:,:,5)); shading flat; title('Non-diatom J term at T = 20°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,5); pcolor(xpar,ChlC,fjld(:,:,1)); shading flat; title('Diatom J term at T = -2°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,6); pcolor(xpar,ChlC,fjld(:,:,2)); shading flat; title('Diatom J term at T = 0°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,7); pcolor(xpar,ChlC,fjld(:,:,3)); shading flat; title('Diatom J term at T = 5°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% % subplot(2,4,8); pcolor(xpar,ChlC,fjld(:,:,5)); shading flat; title('Diatom J term at T = 20°C'); xlabel('Irradiance (W/m2)'); ylabel('Chl/C'); cb = colorbar; ylabel(cb, 'J'); caxis([0 1.9]);
% 
% 
%     % Plot temperature depedence for Jd-irradiance on 2D plots (C/Chl=90 = satellite assumed constant)
%     figure(2);
%     subplot(5,1,1); plot(fjld(:,91,1),'b'); hold on; plot(fjln(:,91,1),'r'); title('J at T=-2°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
%     subplot(5,1,2); plot(fjld(:,91,2),'b'); hold on; plot(fjln(:,91,2),'r'); title('J at T=0°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
%     subplot(5,1,3); plot(fjld(:,91,3),'b'); hold on; plot(fjln(:,91,3),'r'); title('J at T=5°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
%     subplot(5,1,4); plot(fjld(:,91,4),'b'); hold on; plot(fjln(:,91,4),'r'); title('J at T=10°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
%     subplot(5,1,5); plot(fjld(:,91,5),'b'); hold on; plot(fjln(:,91,5),'r'); title('J at T=20°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');

    % Figure 2b - plotting T = -2 and T = 20 curves on same 
    figure(2);
    suptitle('Light-limitation curves for model diatoms and non-diatoms')
    % diatom
    subplot(2,1,1); plot(fjld(:,91,1),'b'); hold on; plot(fjld(:,91,5),'g'); title('Diatom J term ','FontWeight','bold'); legend('Diatom at -2°C', 'Diatom at 20°C'); ylabel('J (dimensionless)'); xlabel('Irradiance (W/m2)'); %xlim([0 10]);
    % non-diatom
    subplot(2,1,2); plot(fjln(:,91,1),'r'); hold on; plot(fjln(:,91,5),'k'); title('Non-diatom J term','FontWeight','bold'); legend('Non-diatom at -2°C', 'Non-diatom at 20°C'); ylabel('J (dimensionless)'); xlabel('Irradiance (W/m2)'); %xlim([0 10]);
%     % both on one plot
%     plot(fjld(:,91,1),'b'); hold on; plot(fjld(:,91,5),'g'); plot(fjln(:,91,1),'r'); plot(fjln(:,91,5),'k'); title('Diatom J term ','FontWeight','bold'); legend('Diatom at -2°C', 'Diatom at 20°C','Non-diatom at -2°C','Non-diatom at 20°C'); ylabel('J'); xlabel('Irradiance (W/m2)'); 

    
    % print -depsc -painters Jterm.eps
    
%     % Plot temperature dependence for Jd-C/Chl on 2D plots (I=50W/m2)
%     figure(3);
%     subplot(4,1,1); plot(ChlC,fjld(51,:,1),'b'); hold on; plot(ChlC,fjln(51,:,1),'r'); title('J at T=-2°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Chl/C');
%     subplot(4,1,2); plot(ChlC,fjld(51,:,2),'b'); hold on; plot(ChlC,fjln(51,:,2),'r'); title('J at T=0°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Chl/C');
%     subplot(4,1,3); plot(ChlC,fjld(51,:,3),'b'); hold on; plot(ChlC,fjln(51,:,3),'r'); title('J at T=5°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Chl/C');
%     subplot(4,1,4); plot(ChlC,fjld(51,:,5),'b'); hold on; plot(ChlC,fjln(51,:,5),'r'); title('J at T=20°C','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Chl/C');
%     
    