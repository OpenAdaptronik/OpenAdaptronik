function [w0, wd]=transff(m, k, d)
    

w0 = (k/m)^.5;
delta = d/(2*m);
if delta >= w0 then wd= 0; return
end
wd = (w0^2 - delta^2)^.5;

w = (wd/2):.01:(wd*3);
s = w*%i;

//figure(3); plot(w,abs(s));

frf = (k + s*d + s.^2*m).^(-1);


f1 = figure('Visible', 'off');

subplot(211);
plot(w, abs(frf));
subplot(212);
plot(w, atan(imag(frf),real(frf)));

xs2png(f1, 'funplot.png');

f=file("open", "eigf.txt", "unknown")
write(f, [w0, wd]);



endfunction
