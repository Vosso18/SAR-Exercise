clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked')

run('minisar_init.m');

figure(1)
    imagesc(abs(raw))
    colormap('jet')
    colorbar()
    xlabel('Azimuth (Samples)')
    ylabel('Range (Samples)')

figure(2)
    imagesc(angle(raw))
    colormap('jet')
    colorbar()
    title('Imaginary')
    xlabel('Azimuth (Samples)')
    ylabel('Range (Samples)')


%% Task 1 (Range and azimuth extent)

% Azimuth resolution
c=3e8;  
lambda=c/sys.fc; 
l=2*sys.azdx;  % Size of azimuth reference cell
roff=c*sys.del/2;  % Range offset
Laz1=lambda*roff/l;  % Azimuth resolution
azpixellenght=sys.azdx;  % Size of azimuth pixel
azpixels=Laz1/azpixellenght  % Number of azimuth pixels

% Range resolution
rangepixellenght=c/(2*sys.fs);  % Size of range pixel
Resnocompr=c/2*sys.pule;  % Range resolution without compression
rangepixels=Resnocompr/rangepixellenght  % Number of range pixels

%% Task 2 (Range compression)


u=@(t) exp(1j*pi*(sys.pubw/sys.pule)*(t).^2);  % Chirp function

N=round(sys.pule*sys.fs*0.5);  
L=2*N+1;  
time=(-N:N)/sys.fs;  

chirp=u(time);  


M=2048;  % Total number of samples in range direction
S=fft(raw);  % Compute FFT of raw data
w=zeros(1,M);  % Weight vector
w(1:L)=(u(time));  
W=fft(w);  
Wconj=transpose(conj(W)); 

ymatrix=zeros(2048,1024);  % Processed data
ymatrixcut=zeros(M-L,1024);  % Invalid samples removed
ymatrixclean=zeros(2048,1024);  % Invalid samples set to zero

for i=1:1024
    y=ifft(S(:,i).*Wconj);
    ymatrix(:,i)=y;

    % Remove invalid samples
    ycut=y(1:M-L);
    ymatrixcut(:,i)=ycut;

    %set invalid samples to zero
    yclean=y;
    yclean(M-L+1:end)=0;
    ymatrixclean(:,i)=yclean;
end

amplcmpr=abs(ymatrix);
phasecmpr=angle(ymatrix);

amplcmprclean=abs(ymatrixclean);
phasecmprclean=angle(ymatrixclean);

figure(4)
imagesc(10*log10(amplcmpr))
colormap("jet")
xlabel('Azimuth (Samples)')
ylabel('Range (Samples)')

figure(5)
imagesc(10*log10(amplcmprclean))
colormap("jet")
xlabel('Azimuth (Samples)')
ylabel('Range (Samples)')


%% Task 3 - Azimuth compression

finalmatrix = zeros(263,1024); % Final matrix


for i = 1:size(ymatrixcut,1) % Over coloumns
    % Chirp signal
    r0 = (i - 1) * rangepixellenght + roff;
    Laz = lambda * r0 * 0.5 / sys.azdx; % Range in m
    K = - 2/(lambda * r0); % Chirp rate
    x = linspace(-Laz/2,Laz/2,Laz); % Length
    u2 = @(x) exp(1j * pi * (-2 / (lambda * r0)) * (x).^2); % Chirp gen
    Naz = round(Laz / (2 * sys.azdx)); % Azimuth range spacing
    L2 = 2 * Naz + 1; % New spacing for start
    space = (-Naz:Naz) * sys.azdx ; % Spacing
    
    % Padding 
    M = size(ymatrixcut,2); % Size of correlation 
    w = zeros(1,M); % Zero pad
    
    
    w(end-Naz+1:end) = (u2(space(1:Naz))); % Begin and end
    w(1:Naz+1) = (u2(space(Naz+1:end)));  % Middle window

    % fft
    W = fft(w); 
    Wconj = conj(W); 
    S = fft(ymatrixcut(i,:)); 
    y = ifft(S .* Wconj);

    finalmatrix(i,:) = y;
end

figure(6)
    imagesc(10*log10(abs(finalmatrix)))
    axis image
    colormap("jet")
    colorbar()

% Chirp signal

figure(7)
    plot(space,real(u2(space)))

% Remove invalid pixels
r0max = (size(ymatrixcut,1)-1) * rangepixellenght + roff;
Lazmax = lambda * r0max * 0.5 / sys.azdx; 
Nazmax = round(Lazmax / (2 * sys.azdx)); 
finalmatrix(:,end-Nazmax:end) = [];
finalmatrix(:,1:Nazmax) = [];


figure(8)
    imagesc(10*log10(abs(finalmatrix)))
    colormap("jet")
    cb = colorbar;
    ylabel(cb, 'Signal Power [dB]');

%% Task 4 - Upsampling

Final_Chop = finalmatrix(121:141, 421:441);
osf = 64;

rows = size(Final_Chop,1);
coloumns = size(Final_Chop,2);

upsamp_row = interpft(Final_Chop, osf*rows,2);      
upsamp_coloumn = interpft(upsamp_row, osf*coloumns,1);      

pow_xcorr = abs(upsamp_coloumn).^2;

figure(9)
    imagesc(pow_xcorr)
    colormap("default")
    colorbar()

%% Task 5 - PSLR and -3 dB

Max = max(pow_xcorr,[],"all");
[ind_x, ind_y] = find(pow_xcorr == Max);

pow_row = pow_xcorr(ind_x,:)/Max;
pow_coloumn = pow_xcorr(:,ind_y)/Max;

pow_row_dB = 10*log10(pow_row); 
pow_col_dB = 10*log10(pow_coloumn);


figure(10)
    plot(pow_row_dB)
    hold on
    yline(-15.882, 'r--', 'PSLR = -15.88 dB', 'LabelHorizontalAlignment','left')
    yline(-3, '--k', '-3 dB = 1.522 m', 'LabelHorizontalAlignment','left')
    hold off
    ylim([-50 0])

figure(11)
    plot(pow_col_dB)
    hold on
    yline(-14.696, 'r--', 'PSLR = -14.70 dB', 'LabelHorizontalAlignment','left')
    yline(-3, '--k', '-3 dB = 1.382 m', 'LabelHorizontalAlignment','left')
    hold off
    ylim([-50 0])


val = -3;

more_pow_row = pow_row_dB > val;
pow_pow_row = find(more_pow_row>0);

Meas_row = length(pow_pow_row) / osf * rangepixellenght


more_pow_col = pow_col_dB > val;
pow_pow_col = find(more_pow_col>0);

Meas_col = length(pow_pow_col) / osf * rangepixellenght


% Theorectical value
th_resolution = (0.886 / sys.pubw * sys.fs) * rangepixellenght


% expect 1.5 and -13.26 (autocorrelation of a sinc function)


