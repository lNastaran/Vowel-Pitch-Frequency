%reading a sound file sampling it by 'Fs' rate into 'original' array
[original, Fs] = audioread('a.wav');
fprintf('Fs : %d\n', Fs);
FrameNum = 123;
%'Fs' sample in each frame,'N' sample in 0.025 seconds
N = 0.025 * Fs;
fprintf('N : %d\n', N);

%'M' --> shift length
M = 0.4 * N;
fprintf('M : %d\n', M);
%number of frames with 'N' length and 'M' shift
NumOfFrames = floor((length(original)-(N-M))/M) + 1;
FramedY = zeros(N, NumOfFrames);
for i = 1:NumOfFrames - 1
    os = ((i - 1) * M) + 1;
    FramedY(:,i) = original(os : os + N - 1);
end

%the last frame is not full and we fill it with zeros
FramedY(1:(length(original) - (i * M)),i+1) = original(((i * M) + 1):length(original));

%windowing
%hamming works better
window = hamming(N);

%windowing each frame by multiply each voice signal value with corresponding window value
WindowedY = zeros(N, NumOfFrames);
for k = 1:NumOfFrames
    WindowedY(:,k) = FramedY(:,k) .* window;
end

%Preemphasis
alpha = 0.9;
for k = 1:NumOfFrames
    for j = 2:N
        WindowedY(j,k) = WindowedY(j,k) - (alpha * WindowedY(j - 1,k));
    end
end

%which frame is vowel
%after that we will find pitch freq on a random frame
for k = 1:NumOfFrames
    DC = 0;
    %calculating DC value by finde the mean of signal
    
    for c = 1:N
            DC = DC + WindowedY(c,k);
    end
    DC = DC / N;
   
    %set DC to zero
    ZeroOffsetWindowedY(:,k) = WindowedY(:,k) - DC;
end

%Energy of each frame
temp = ZeroOffsetWindowedY .* ZeroOffsetWindowedY;
energy = sum(temp)/N;

%range variable specifies how many frame of start and end of filewe assum zero
range = 4;
ESilence = (sum(energy(1:range))+sum(energy(length(energy) - range + 1:length(energy))))/(range*2);

%applying a threshhold to remove silence part of signal
SilenceThresh = 2.5 * ESilence;



%applying fast fourie transform on windowed signal
WindowedYFreq = fft(WindowedY);

%logarithm of absolout of fast fourie of windowed signal
Logabs = log(abs(WindowedYFreq));

%reverse fast foure on the last value to achive cepstrum
%cepstrum : fft() --> Log(abs()) --> ifft()
%at last we need real part of it
cepstrum=real(ifft(Logabs));

%in order to get pitch frequency we need to get high time lifter from cepstrum
High_time_cepstrum = zeros(size(cepstrum(:,1),1)/2, NumOfFrames);
%threshold=20
for k1 = 1:NumOfFrames
    mul = 0;
    for k2 = 1:(size(cepstrum(:,k1),1)/2)
        if (k2 == 20)
            mul = 1;
        end
        High_time_cepstrum(k2,k1) = cepstrum(k2,k1) .* mul;
    end
end


%peaks location
[pks, locs]= findpeaks(High_time_cepstrum(:,FrameNum),"DoubleSided");
[I1,I2] = max(pks);

x = [1:200];
x = Fs./x;

%sample to hz
pitch = Fs./locs(I2);

subplot(5, 1, 1);
plot(original);
TitleStr = sprintf('Original Audio');
title(TitleStr);

subplot(5, 1, 2);
plot(WindowedY(:,FrameNum));
TitleStr = sprintf('Frame number %d after windowing', FrameNum);
title(TitleStr);

subplot(5, 1, 3);
plot(abs(WindowedYFreq(:,FrameNum)));
TitleStr = sprintf('Fast Fourie Transform of frame number %d', FrameNum);
title(TitleStr);

subplot(5, 1, 4);
plot(cepstrum(2:400,FrameNum));
TitleStr = sprintf('Cepstrum of frame number %d', FrameNum);
title(TitleStr);

subplot(5, 1, 5);
plot(High_time_cepstrum(1:200,FrameNum));
TitleStr = sprintf('High Time Cepstrum of frame number %d', FrameNum);
title(TitleStr);

fprintf("Pitch frequency = %f Hz\n",pitch);

set(gcf, 'Position', get(0, 'Screensize'));