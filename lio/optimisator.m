global errs;

errs = [];

x0 = zeros(49*49,1);

x = fminsearch(@fmin, x0);



function s = fmin(x)
    global errs;
    
    s = sum(sum(takeImage(x)));
    
    errs = [errs s];
end