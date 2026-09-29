function rgb=hex2rgb(value)
value=char(value);if startsWith(value,'#'),value=value(2:end);end
assert(numel(value)==6,'Expected a six-digit hexadecimal color.');
rgb=reshape(sscanf(value,'%2x'),1,3)/255;
end
