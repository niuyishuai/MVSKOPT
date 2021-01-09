function fval = func_subs(f,var,val,modtype)
    % subs var of f as val for model of type modtype
    switch modtype
        case "yalmip"
            fval = replace(f,var,val);
        case "polylab"
            fval = f.eval(val);
        case "sym"
            fval = double(subs(f,var,val));
        case "optimvar"
            x.x=val;
            fval = evaluate(f,x);
        case "sostools"
            fval = double(subs(f,var,val));
        otherwise
             fval = double(subs(f,var,val));
    end
end