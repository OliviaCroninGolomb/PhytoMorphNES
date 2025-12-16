function [wvln,ap,aph,ad,abs_ap,abs_ad] = apspecproc_Dec2012(apfilename, adfilename, volfilt, filterarea, zerowvlnmin, zerowvlnmax)
    % apspecproc_Dec2012: Process raw absorbance data into absorption coefficients
    %
    % Inputs:
    %   apfilename  - Path to ap file (particle absorption).
    %   adfilename  - Path to ad file (detritus absorption).
    %   volfilt     - Volume filtered in cubic meters.
    %   filterarea  - Area of filter in square meters.
    %   zerowvlnmin - Minimum wavelength for zero offset.
    %   zerowvlnmax - Maximum wavelength for zero offset.
    %
    % Outputs:
    %   wvln        - Wavelengths.
    %   ap          - Particle absorption coefficients.
    %   aph         - Phytoplankton absorption coefficients.
    %   ad          - Detritus absorption coefficients.
    %   abs_ap      - Raw absorbance for ap.
    %   abs_ad      - Raw absorbance for ad.

    % Constants for Mitchell et al. (1990)
    Mitcha = 0.392;
    Mitchb = 0.655;

    % Read ap data (replace with new format)
    ap_data = readtable(apfilename, 'FileType', 'text', 'HeaderLines', 86);
    wvln = ap_data{:, 1};          % Wavelengths
    abs_ap = ap_data{:, 2};        % Raw absorbance
    
    % Zero offset correction
    zero_range = (wvln >= zerowvlnmin & wvln <= zerowvlnmax);
    ap_zero = mean(abs_ap(zero_range));
    ap_corrected = abs_ap - ap_zero;
    
    % Compute optical density
    ap_od = Mitcha * ap_corrected + Mitchb * ap_corrected.^2;
    
    % Compute absorption coefficients
    ap = 2.3 * ap_od / (volfilt / filterarea);

    % Read ad data if provided
    if ~isempty(adfilename)
        ad_data = readtable(adfilename, 'FileType', 'text', 'HeaderLines', 86);
        abs_ad = ad_data{:, 2};
        ad_zero = mean(abs_ad(zero_range));
        ad_corrected = abs_ad - ad_zero;
        ad_od = Mitcha * ad_corrected + Mitchb * ad_corrected.^2;
        ad = 2.3 * ad_od / (volfilt / filterarea);
    else
        abs_ad = [];
        ad = zeros(size(ap)) - 9999;  % Placeholder for no ad data
    end

    % Phytoplankton absorption
    aph = ap - ad;
end
