function switch_template(template)

switch template
    case 'cCARLIN'
        disp("----------Use cCARLIN template---------")
        copyfile('@CARLIN_def/CARLIN_def_cCARLIN.m','@CARLIN_def/CARLIN_def.m')
    case 'Tigre'
        disp("----------Use Tigre template----------")
        copyfile('@CARLIN_def/CARLIN_def_Tigre.m','@CARLIN_def/CARLIN_def.m')
    case 'Tigre_2022'
        disp("----------Use Tigre template----------")
        copyfile('@CARLIN_def/CARLIN_def_Tigre_2022.m','@CARLIN_def/CARLIN_def.m')
    case 'Tigre_2022_v2'
        disp("----------Use Tigre template----------")
        copyfile('@CARLIN_def/CARLIN_def_Tigre_2022_v2.m','@CARLIN_def/CARLIN_def.m')
    case 'Rosa'
        disp("----------Use Rosa template----------")
        copyfile('@CARLIN_def/CARLIN_def_Rosa.m','@CARLIN_def/CARLIN_def.m')
    case 'Rosa_v2'
        disp("----------Use Rosa template----------")
        copyfile('@CARLIN_def/CARLIN_def_Rosa_v2.m','@CARLIN_def/CARLIN_def.m')
        
    otherwise
        error('Input should be Tigre, Tigre_2022, Tigre_2022_v2, cCARLIN, Rosa, or Rosa_v2')
end


