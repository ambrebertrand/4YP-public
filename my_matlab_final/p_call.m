% Pushbutton callback
function p_call(src, event, h)
    vals = get(h.c,'Value');
    checked = find([vals{:}]);
    if isempty(checked)
        checked = 'none';
    end
    assignin('base','checked',checked)
    close(h.f) 
end
