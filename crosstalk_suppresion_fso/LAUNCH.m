% Launcher

poolnum = 1;

params = params_setup(0, 0);
SNR = params.SNR_dB_vec;
r0 = params.r0_vec;

MAX_i = length(SNR);

for i = 1:MAX_i
    Parallel_Experiment_conduct (SNR(i), r0(i), poolnum)
end

poolobj = gcp('nocreate');
delete(poolobj);