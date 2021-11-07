%reading a sound file sampling it by 'Fs' rate into 'original' array
[original, Fs] = audioread('a.wav');
FrameNum = 123;
%'Fs' sample in each frame,'N' sample in 0.025 seconds
N = 0.025 * Fs;

%'M' --> shift length
M = 0.4 * N;
%number of frames with 'N' length and 'M' shift
NumOfFrames = floor((length(original)-(N-M))/M) + 1;
FramedY = zeros(N, NumOfFrames);
for i = 1:NumOfFrames - 1
    os = ((i - 1) * M) + 1;
    FramedY(:,i) = original(os : os + N - 1);
end

%the last frame is not full and we fill it with zeros
FramedY(1:(length(original) - (i * M)),i+1) = original(((i * M) + 1):length(original));

%choose windowing
window=0;%hamming works better
if window == 1
    wind = rectwin(N);
else
    wind = hamming(N);
end

%windowing each frame by multiply each voice signal value with corresponding window value
WindowedY = zeros(N, NumOfFrames);
for k = 1:NumOfFrames
    WindowedY(:,k) = FramedY(:,k) .* wind;
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

%in order to get formants we need to get low time lifter from cepstrum
Low_time_cepstrum = zeros(size(cepstrum(:,1),1)/2, NumOfFrames);
%threshold=20
for k1 = 1:NumOfFrames
    mul = 1;
    for k2 = 1:(size(cepstrum(:,k1),1)/2)
        if (k2 == 15)
            mul = 0;
        end
        Low_time_cepstrum(k2,k1) = cepstrum(k2,k1) .* mul;
    end
end


Ch = fft(Low_time_cepstrum(:,FrameNum));
Ch = log(abs(Ch));
OrderedCh=sort(Ch);
[peak,location]=max(Ch);
peak2=OrderedCh(2);
peak3=OrderedCh(3);
location2=find(Ch==peak2)(1);
location3=find(Ch==peak3)(1);

x = [1:200];
x = Fs./x;

%sample to hz
formant1 = Fs./location;
formant2 = Fs./location2;
formant3 = Fs./location3;

subplot(5, 1, 1);
plot(original);
TitleStr = sprintf('Original Audio Time Signal');
title(TitleStr);

subplot(5, 1, 2);
plot(WindowedY(:,FrameNum));
TitleStr = sprintf('Frame number %d in Hamming', FrameNum);
title(TitleStr);

subplot(5, 1, 3);
plot(abs(WindowedYFreq(:,FrameNum)));
TitleStr = sprintf('Fast Fourie Transform of %d frame', FrameNum);
title(TitleStr);

subplot(5, 1, 4);
plot(cepstrum(2:400,FrameNum));
TitleStr = sprintf('Log absoulote of %d ftame fft or Cepstrom', FrameNum);
title(TitleStr);

subplot(5, 1, 5);
plot(Low_time_cepstrum(1:200,FrameNum));
TitleStr = sprintf('Low Time Cepstrom of frame %d', FrameNum);
title(TitleStr);

fprintf("first formant = %f hz\n",formant1);
fprintf("second formant = %f hz\n",formant2);
fprintf("third formant = %f hz\n",formant3);


set(gcf, 'Position', get(0, 'Screensize'));