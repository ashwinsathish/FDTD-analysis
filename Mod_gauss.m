fc = 7.5e9; % Center frequency
bw = 9e9;
fm = bw/2;
sigma = bw/(2*sqrt(2*log(2)));

%************** Iteration loop ****************

   
for n=1:100         
    time=(n-1)*dt;  
        
%------- Launching the signal ---------

      m = cos(2*pi*fm*time);

      % Gaussian envelope
      g = exp(-(time-t0).^2/(2*sigma^2));

      % modulated Gaussian pulse
      pulse = g.*m;
end

figure;
plot(time, pulse);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modulated Gaussian Pulse with Sine Wave Modulation');
grid on;
