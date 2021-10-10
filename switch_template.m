function switch_template(template)

switch template
    case 'cCARLIN'
        copyfile('@CARLIN_def/CARLIN_def_cCARLIN.m','@CARLIN_def/CARLIN_def.m')
    case 'Tigre'
        copyfile('@CARLIN_def/CARLIN_def_Tigre.m','@CARLIN_def/CARLIN_def.m')
    otherwise
        error('Input should be Tigre or cCARLIN')
end


