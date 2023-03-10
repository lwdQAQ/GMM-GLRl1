function delta_t = calc_time_step_adaptive(eval_obj_fun, update_params, ...
    param, params, der_param, delta_t, delta_t0, Y)
% If the attempted delta_t is less than the lower limit, reset it to the
% lower limit.
if delta_t < delta_t0
    delta_t = delta_t0;
end

val_ori = eval_obj_fun(param, Y, params);
param_new = update_params(param, der_param, delta_t);
val_new = eval_obj_fun(param_new, Y, params);

if val_new < val_ori
    val_old = val_new;
    while 1
        delta_t = delta_t * 10;
        param_new = update_params(param, der_param, delta_t);
        val = eval_obj_fun(param_new, Y, params);
        if val < val_old
            val_old = val;
        else
            delta_t = delta_t / 10;
            break;
        end
    end
else
    while delta_t >= delta_t0
        delta_t = delta_t / 10;
        param_new = update_params(param, der_param, delta_t);
        val = eval_obj_fun(param_new, Y, params);
        if val < val_ori
            break;
        end
    end
    if delta_t < delta_t0
        delta_t = 0;
    end
end


