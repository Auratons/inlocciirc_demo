function env = environment()
    if contains(pwd, 'datagrid')
        env = 'cmp';
    elseif contains(pwd, 'repos')
        env = 'laptop';
    else
        env = 'ciirc';
    end
end