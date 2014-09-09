% J sensitivty to alpha (P-I slope scaling factor that feeds into J)

% EDIT: Script not tested since edited for ChlC error spotted in
% plottingJ.m (21/03/2014)


% J(I,CChl,alpha)

alpn = [2 15 20 30 50];
alpd = [2 11.25 20 30 50];

for i = 1:length(alpn)
    
% Constants

loc_T = double(0); %double(-2); % water temperature
% loc_T = Temp(i); % For testing J sensitivity to T
xvpn = double(0.53); % xvpn = Vpn in paper - max nd phytoplankton growth rate at 0°C. = 0.53 (Yool2010 table 1)
xvpd = double(0.50); % xvpd = Vpd in paper - max d phytoplankton growth rate at 0°C. = 0.50 (Yool2010 table 1)


% Temperature depedence

fun_T = 1.066 .* exp(loc_T); % Yool2010 eq14
xvpnT = xvpn .* fun_T; 
xvpdT = xvpd .* fun_T; 



% Light field and photodependence (Chl/biomass)

% alphan = double(15.0); alphad = double(11.25); % Non-diatom and diatom chl-specific initial slope of P-I curve (Yool2010 Table1)
alphan = alpn(i); alphad = alpd(i);

% edited section
    ChlC   = 0:(0.05/100):0.05; % 90; % biomass/Chl: in situ published range = 24-1250 (OASIS_MATRAI_2013 slides), model range = 50-200 (model enforced min of 20, Yool 2010 Table 1), satellite = 90.
    xpar   = 0:1:100; % Matsuoka2009 fig2a = 0-359 W/m2 (sea surface).

    faln = alphan * CChl; % Yool2010 eq13
    fald = alphad * CChl; % Yool2010 eq20

    faln (faln > 0.05) = 0.05; % temporary! - it assumes Yool2013 table 1 error (that faln(max) = 0.05 (not ChlC(max)).
    fald (fald > 0.05) = 0.05; % temporary! - it assumes Yool2013 table 1 error


pos = 1;
    for j = 1:length(faln);
    
        % non-diatoms
        fchn1 = (xvpnT.^2) + ((faln(j).^2) .* (xpar.^2)); % Yool2010 eq15 (denominator)
        fchn  = xvpnT ./ sqrt(fchn1);
        fjln(:,pos,i)  = fchn .* faln(j) .* xpar; % Yool2010 eq 15 (full)

        % diatoms
        fchd1 = (xvpdT.^2) + ((fald(j).^2) .* (xpar.^2)); % Yool2010 eq22 (denominator)
        fchd  = xvpdT ./ sqrt(fchd1);
        fjld(:,pos,i)  = fchd .* fald(j) .* xpar; % Yool2010 eq22 (full)

        pos = pos + 1;
    end
    
end


% Subplot grid of temperature dependence on J term 3D plot
figure(1);
subplot(2,5,1); pcolor(fjln(:,:,1)); shading flat; title('Non-diatom J term at alpha = 2'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,2); pcolor(fjln(:,:,2)); shading flat; title('Non-diatom J term at alpha = 15 (Yool)'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,3); pcolor(fjln(:,:,3)); shading flat; title('Non-diatom J term at alpha = 20'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,4); pcolor(fjln(:,:,4)); shading flat; title('Non-diatom J term at alpha = 30'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,5); pcolor(fjln(:,:,5)); shading flat; title('Non-diatom J term at alpha = 50'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,6); pcolor(fjld(:,:,1)); shading flat; title('Diatom J term at alpha = 2'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,7); pcolor(fjld(:,:,2)); shading flat; title('Diatom J term at alpha = 11.25 (Yool)'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,8); pcolor(fjld(:,:,3)); shading flat; title('Diatom J term at alpha = 20'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,9); pcolor(fjld(:,:,4)); shading flat; title('Diatom J term at alpha = 30'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');
subplot(2,5,10); pcolor(fjld(:,:,5)); shading flat; title('Diatom J term at alpha = 50'); xlabel('Irradiance (W/m2)'); ylabel('C/Chl'); cb = colorbar; ylabel(cb, 'J');


    % Plot temperature depedence for Jd-irradiance on 2D plots (C/Chl=90)
    figure(3);
    subplot(5,1,1); plot(fjld(:,91,1),'b'); hold on; plot(fjln(:,91,1),'r'); title('J at alpha = 2','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
    subplot(5,1,2); plot(fjld(:,91,2),'b'); hold on; plot(fjln(:,91,2),'r'); title('J at alphan = 15, alphad = 11.25','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
    subplot(5,1,3); plot(fjld(:,91,3),'b'); hold on; plot(fjln(:,91,3),'r'); title('J at alpha = 20','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
    subplot(5,1,4); plot(fjld(:,91,4),'b'); hold on; plot(fjln(:,91,4),'r'); title('J at alpha = 30','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');
    subplot(5,1,5); plot(fjld(:,91,5),'b'); hold on; plot(fjln(:,91,5),'r'); title('J at alpha = 50','FontWeight','bold'); legend('Diatom J','Non-diatom J'); ylabel('J'); xlabel('Irradiance (W/m2)');