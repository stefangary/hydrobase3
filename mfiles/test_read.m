% script to test Hydrobase-M read functions ...

fidin = hb_open_file('', '5603', '.btl');

err = 0;
nsta = 0;
while (err ~= -1)
    nsta = nsta + 1;
    [profiles(nsta),err] = hb_get_profile(fidin);
end

fclose(fidin);
clear fidin;



